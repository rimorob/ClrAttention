# Fig. 1 (method): CLR as attention. Big-picture strip (compendium -> four steps -> module map) above a
# worked toy example that traces one gene pair (i, j) and a third gene k through the four steps.
# Double-column figure ONLY (7.2 in wide); it is not designed for a single column.
#
# Layout (rev. 3): strip (input, steps with 'published CLR 2007 / added here' bracket, output) over a
# 2 x 2 grid of numbered panels read like text (1, 2 / 3, 4). Panel 3 starts again from S, labelled
# 'S, from step 2' (no long connector: the numbered chips carry the order). Layout inspiration: the
# PaperBanana draft of 30 Oct (layout only; its matrices and numbers were invented, everything here is
# computed).
#
# REAL TOY: every matrix is computed with the package's own functions on a small synthetic
# expression matrix with planted structure (a chain module, a one-factor module, noise genes):
#   D     bspline_mi() via ClrAttention$estimate_mi()        (diagonal zeroed as clr_calibrate does)
#   CLR   clr_calibrate(method = "normal", combine = "stouffer") via $calibrate()
#   A     .softmax_rows(S, keff, cap) via $build_operator(softmax_keff, softmax_cap)
#   P^t   P = (1-alpha) I + alpha A, repeated product, K_t = (P^t + P^t')/2: the lines of
#         analysis/baselines.R (the walk as evaluated in the paper); ASSERTED below by evaluating those
#         very lines of baselines.R on the toy and by an independent repeated-squaring computation.
# The z-score curves are recomputed from the same formulas as src/core/clr_core.cpp and ASSERTED equal
# to the clr_calibrate() output. Diagonal convention (clr_core.cpp, following clr.m): the diagonal of D
# is set to 0 and KEPT among the G values that give mu_i and sigma_i (sample sd, n - 1); it is zeroed
# again in S, and drawn grey in the matrices.
# Toy parameters differ from the real method only where a 14-gene matrix forces it:
#   cap 5 (real 50), k_eff 3 (real 10; with <= k_eff candidates a row would be uniform), t = 1/3/9 shown.
#
# Colour key as in analysis/plot_theme.R: blue = dependence (MI), orange = CLR normalization,
# violet = softmax, magenta = depth / walk, grey = neutral. Palette decision (review F12): gene identity
# no longer uses key colours. Gene i = solid black outline, gene j = dashed black outline, both on a
# grey band; k = triangle. So blue stays 'dependence' and magenta stays 'walk' everywhere.
#
# Output: docs/whitepapers/figures/fig1_method.{png,pdf}
# Run:    nice -n 19 Rscript analysis/fig_method.R

suppressPackageStartupMessages({ library(grid); library(clr); library(Matrix) })
source("analysis/plot_theme.R")
out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"

# ---- toy data and the real pipeline --------------------------------------------------------------
SEED <- 7; G <- 14; N <- 60
CAP <- 5L; KEFF <- 3; ALPHA <- 0.5; TS <- c(1L, 3L, 9L); TUSED <- 9L
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

# ---- assertions: the picture cannot drift from the code -------------------------------------------
zraw <- c(zi$z[J], zj$z[I])
stopifnot(abs((pmax(zraw[1], 0) + pmax(zraw[2], 0)) / sqrt(2) - S[I, J]) < 1e-10)
Zfull <- t(apply(D, 1, function(r) (r - mean(r)) / sd(r))); Zfull[Zfull < 0] <- 0
stopifnot(max(abs((Zfull + t(Zfull)) / sqrt(2) - replace(S, cbind(1:G, 1:G), 0))) < 1e-10)
# printed '(2.2 + 1.8)/sqrt(2) = 2.8' is arithmetically consistent with the rounded z values
stopifnot(round((round(zi$z[J], 1) + round(zj$z[I], 1)) / sqrt(2), 1) == round(S[I, J], 1),
          zi$z[J] > 0, zj$z[I] > 0)                                 # no clipping for this pair
# diagonal convention: zeroed diagonal enters mu and sigma (G values, not G - 1)
stopifnot(length(D[I, ]) == G, D[I, I] == 0,
          abs(zi$mu - mean(D[I, -I]) * (G - 1) / G) < 1e-12)
# softmax: row sums, k_eff for ALL rows with more than k_eff candidates, cap binds for row i
cand <- rowSums(A > 0)
ent <- function(a) exp(-sum(a[a > 0] * log(a[a > 0])))
rows_ok <- which(cand > KEFF)
stopifnot(all(abs(rowSums(A) - 1) < 1e-9), length(rows_ok) == G,
          all(abs(vapply(rows_ok, function(r) ent(A[r, ]), 0) - KEFF) < 1e-4),
          cand[I] == CAP, all(cand <= CAP))
# gene k: positive CLR score but the 6th positive entry of row i, cut by the toy cap -> A_ik = 0 (and A_ki = 0)
pos_i <- sort(S[I, S[I, ] > 0], decreasing = TRUE)
stopifnot(S[I, K] > 0, A[I, K] == 0, A[K, I] == 0, length(pos_i) > CAP,
          match(S[I, K], pos_i) == CAP + 1L, A[I, J] == max(A[I, ]))
# walk: depth ordering at (i, k) and K_3 < K_9, K_1 = 0 exactly
stopifnot(Kt[["1"]][I, K] == 0, Kt[["1"]][I, K] < Kt[["3"]][I, K], Kt[["3"]][I, K] > 0,
          Kt[["3"]][I, K] < Kt[["9"]][I, K])
# walk equals (a) the lines of analysis/baselines.R evaluated on this toy and (b) an independent
# repeated-squaring computation of P^t
bl <- readLines("analysis/baselines.R")
i0 <- grep("^P <- Matrix::Diagonal\\(G, 1 - opt\\$alpha\\)", bl)
i1 <- grep("^St <- \\(Pt \\+ t\\(Pt\\)\\) / 2", bl)
stopifnot(length(i0) == 1, length(i1) == 1, i1 > i0, i1 - i0 < 6)
env <- new.env(); env$G <- G; env$opt <- list(alpha = ALPHA, t = TUSED); env$f <- fit
eval(parse(text = bl[i0:i1]), env)
stopifnot(max(abs(env$St - Kt[[as.character(TUSED)]])) < 1e-12)
mpow <- function(M, n) { R <- diag(nrow(M)); B <- M; while (n > 0) { if (n %% 2) R <- R %*% B; B <- B %*% B; n <- n %/% 2 }; R }
Pd <- as.matrix(P)
for (tt in TS) { Pt <- mpow(Pd, tt); stopifnot(max(abs((Pt + t(Pt)) / 2 - Kt[[as.character(tt)]])) < 1e-12) }
stopifnot(isSymmetric(Kt[["9"]]), all(abs(rowSums(mpow(Pd, 9)) - 1) < 1e-9))

# ---- 'what the walk buys' (illustrative): co-membership AP, plain CLR versus K_t ---------------
# One seed, 14 genes, PLANTED modules: a sanity illustration of the mechanism, not evidence.
modlab <- integer(G); modlab[chain] <- 1L; modlab[modB] <- 2L
ut <- upper.tri(S)
truth <- (outer(modlab, modlab, function(a, b) a > 0 & a == b))[ut]
ap_one <- function(score, truth) { o <- order(-score, runif(length(score))); t <- truth[o]
  mean(cumsum(t)[t] / which(t)) }
ap <- function(M, reps = 400L) { set.seed(1); sc <- M[ut]               # random tie-breaking, seeded
  mean(replicate(reps, ap_one(sc, truth))) }
stopifnot(abs(ap_one(c(3, 2, 1), c(TRUE, FALSE, TRUE)) - (1 + 2 / 3) / 2) < 1e-12,   # AP unit test
          sum(truth) == choose(length(chain), 2) + choose(length(modB), 2),
          length(unique(Kt[["9"]][ut])) == sum(ut))                      # K_9 has no ties
AP_CLR <- ap(S); AP_K <- vapply(Kt, ap, 0); AP_BASE <- mean(truth)
# NOTE: no assertion that the walk wins; the AP values are printed and drawn as they come out.
cat(sprintf("toy: i=%d j=%d k=%d  D_ij=%.3f  z_ij=%.2f z_ji=%.2f  S_ij=%.2f  A_ij=%.2f  A_ik=%g  K_t[i,k]=%s\n",
            I, J, K, D[I, J], zi$z[J], zj$z[I], S[I, J], A[I, J], A[I, K],
            paste(sprintf("t%d:%.3f", TS, sapply(Kt, function(m) m[I, K])), collapse = " ")))
cat(sprintf("mu_i=%.4f sd_i=%.4f mu_j=%.4f sd_j=%.4f; S_ik=%.3f (rank %d of positives); AP CLR=%.3f K_t=%s base=%.3f\n",
            zi$mu, zi$sd, zj$mu, zj$sd, S[I, K], match(S[I, K], pos_i), AP_CLR,
            paste(sprintf("%.3f", AP_K), collapse = "/"), AP_BASE))
cat("all assertions passed\n")

# ---- drawing helpers (canvas in inches; y measured UP, Y(d) converts a distance from the top) -------
W <- 7.2; H <- 8.0
Y <- function(d) H - d
XD <- 3.5                                          # divider between the two columns of the grid
fill_grey <- "#52514e"
ramp_of <- function(hi, mid) colorRampPalette(c("#fcfcfb", mid, hi))(100)
ramps <- list(D = ramp_of("#184f95", "#9ec5f4"), S = ramp_of("#a8390f", "#f6b595"),
              A = ramp_of("#3a2aa0", "#c9c2f4"), K = ramp_of("#a8345f", "#f5bfd3"),
              X = colorRampPalette(c("#fcfcfb", "#b9b8b3", "#3d3c3a"))(100))
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
rect <- function(x, y, w, h, fill = NA, col = NA, lwd = 1, alpha = 1, lty = 1, just = c("left", "bottom"))
  grid.rect(x, y, w, h, default.units = "native", just = just,
            gp = gpar(fill = fill, col = col, lwd = lwd, alpha = alpha, lty = lty, linejoin = "mitre"))
head_arrow <- function(...) seg(..., arrow = arrow(length = unit(0.075, "inches"), type = "closed"))
big_arrow <- function(x0, x1, y, col, h = 0.17, hl = 0.17) {   # filled block arrow
  grid.polygon(c(x0, x1 - hl, x1 - hl, x1, x1 - hl, x1 - hl, x0),
               c(y + h / 3, y + h / 3, y + h / 1.4, y, y - h / 1.4, y - h / 3, y - h / 3),
               default.units = "native", gp = gpar(fill = col, col = NA))
}
heat <- function(x0, ytop, cell, M, ramp, vmax = max(M), gamma = 0.7) {   # matrix, top-left (x0, ytop)
  n <- nrow(M); r <- rep(seq_len(n), times = n); cc <- rep(seq_len(n), each = n)
  v <- pmin(pmax(as.vector(M) / vmax, 0), 1)^gamma
  fl <- ramp[pmax(1, ceiling(v * 100))]
  fl[r == cc] <- "#ececea"
  rect(x0 + (cc - 1) * cell, ytop - r * cell, cell, cell, fill = fl, col = "white", lwd = 0.25)
}
# identity marks: gene i = solid black frame on a grey band, gene j = dashed black frame on a grey band
bands <- function(x0, ytop, cell, n, i, j, ext = 0.1, i_on = TRUE, j_on = TRUE) {
  if (i_on) { rect(x0 - ext, ytop - i * cell, n * cell + 2 * ext, cell, fill = fill_grey, alpha = 0.2)
    rect(x0, ytop - i * cell, n * cell, cell, col = "black", lwd = 1) }
  if (j_on) { rect(x0 + (j - 1) * cell, ytop - n * cell - ext, cell, n * cell + 2 * ext, fill = fill_grey, alpha = 0.2)
    rect(x0 + (j - 1) * cell, ytop - n * cell, cell, n * cell, col = "black", lwd = 1, lty = 2) }
}
mark <- function(x0, ytop, cell, i, j, col = "black", lwd = 1.7)   # traced cell
  rect(x0 + (j - 1) * cell, ytop - i * cell, cell, cell, col = col, lwd = lwd)
tri <- function(x, y, size = 4.5) grid.points(x, y, default.units = "native", pch = 17,
                                              size = unit(size, "pt"), gp = gpar(col = "black"))
chip <- function(num, x, y, r = 0.105) {
  grid.circle(x, y, r, default.units = "native", gp = gpar(fill = ink[["primary"]], col = NA))
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
panel_head <- function(num, title, x, y, tagtxt, col, right) {
  chip(num, x + 0.105, y)
  txt(title, x + 0.27, y, 9.5, face = "bold", just = "left")
  tag(tagtxt, x, y - 0.27, col)
  for (q in seq_along(right)) txt(right[q], x + 3.35, y + (length(right) - 1) * 0.06 - (q - 1) * 0.12, 7,
                                  ink[["secondary"]], "italic", just = "right")
}
pill <- function(label, x, y, h = 0.2, size = 7, fill = "#f1f0ed", col = ink[["primary"]], face = "plain",
                 border = NA) {
  w <- convertWidth(grobWidth(textGrob(label, gp = gpar(fontfamily = FONT, fontsize = size, fontface = face))),
                    "native", TRUE) + 0.2
  grid.roundrect(x, y, w, h, default.units = "native", r = unit(h / 2, "inches"),
                 gp = gpar(fill = fill, col = border, lwd = 0.8))
  txt(label, x, y, size, col, face)
}

# distribution of one row's MI values on a SHARED MI axis: the normal fit CLR assumes (it extends below
# MI = 0, which is impossible for MI; shaded), mu and sigma marked, z * sigma bracketed
LO <- -0.09; HI <- 0.37
stopifnot(max(D[c(I, J), ]) < HI)
dist_panel <- function(x0, x1, ybase, h, zz, obs_idx, who) {
  vals <- D[who, ]; mu <- zz$mu; s <- zz$sd; obs <- vals[obs_idx]
  xs <- seq(LO, HI, length.out = 300); dens <- dnorm(xs, mu, s); ymax <- dnorm(mu, mu, s)
  X_ <- function(v) x0 + (v - LO) / (HI - LO) * (x1 - x0)
  Y_ <- function(d) ybase + d / ymax * h
  rect(X_(LO), ybase, X_(0) - X_(LO), h * 1.12, fill = ink[["grid"]], alpha = 0.8)        # MI < 0 zone
  grid.polygon(c(X_(xs), X_(HI), X_(LO)), c(Y_(dens), ybase, ybase), default.units = "native",
               gp = gpar(fill = "#9ec5f4", col = NA, alpha = 0.55))
  tl <- xs >= obs
  grid.polygon(c(X_(xs[tl]), X_(HI), X_(obs)), c(Y_(dens[tl]), ybase, ybase),
               default.units = "native", gp = gpar(fill = "#184f95", col = NA, alpha = 0.9))
  grid.lines(X_(xs), Y_(dens), default.units = "native", gp = gpar(col = "#184f95", lwd = 1.2))
  seg(x0, ybase, x1, ybase, ink[["secondary"]], 0.8)
  for (tk in c(0, 0.1, 0.2, 0.3)) { seg(X_(tk), ybase, X_(tk), ybase - 0.04, ink[["secondary"]], 0.7)
    txt(sub("^0\\.", ".", sprintf("%g", tk)), X_(tk), ybase - 0.12, 7, ink[["secondary"]]) }
  others <- setdiff(seq_along(vals), c(who, obs_idx))    # the row's other MI values as a rug
  seg(X_(vals[others]), ybase, X_(vals[others]), ybase + 0.045, "#184f95", 0.6)
  seg(X_(vals[who]), ybase, X_(vals[who]), ybase + 0.045, ink[["muted"]], 1.1)   # zeroed diagonal (grey)
  seg(X_(mu), ybase, X_(mu), Y_(ymax) + 0.03, ink[["secondary"]], 1.1)         # mu (solid grey: dashed is gene j)
  seg(X_(mu), Y_(ymax * exp(-0.5)), X_(mu + s), Y_(ymax * exp(-0.5)), ink[["primary"]], 1.4)   # sigma
  yb <- ybase + h * 1.5                                                            # z * sigma bracket
  seg(X_(mu), yb, X_(obs), yb, ink[["primary"]], 0.9)
  seg(X_(mu), yb - 0.03, X_(mu), yb + 0.03, ink[["primary"]], 0.9)
  seg(X_(obs), yb - 0.03, X_(obs), yb + 0.03, ink[["primary"]], 0.9)
  seg(X_(obs), yb, X_(obs), ybase, ink[["primary"]], 0.6, lty = 3)
  grid.points(X_(obs), ybase, default.units = "native", pch = 19, size = unit(5, "pt"),
              gp = gpar(col = "#184f95"))
  list(X = X_, yb = yb, obs = obs, mu = mu, s = s, Y = Y_, ymax = ymax)
}

# ---- the figure ---------------------------------------------------------------------------------
ORT <- "#b5461a"                                   # darker orange for TEXT only (contrast on white)
draw_fig <- function() {
  grid.newpage(); grid.rect(gp = gpar(fill = "white", col = NA))
  pushViewport(viewport(xscale = c(0, W), yscale = c(0, H), clip = "on"))

  # ================= headline + badge ==============================================================
  txt("CLR as attention: compendium in, gene-association map out", 0.1, Y(0.15), 8, face = "bold", just = "left")
  pill("WORKED EXAMPLE (TOY): not results", W - 1.1, Y(0.15), h = 0.2, size = 7, fill = ink[["primary"]],
       col = "white", face = "bold")

  # ================= big-picture strip ============================================================
  s_top <- 0.33; s_bot <- 1.4
  grid.roundrect(0.1, Y(s_bot), W - 0.2, s_bot - s_top, default.units = "native", just = c("left", "bottom"),
                 r = unit(0.05, "inches"), gp = gpar(fill = "#f7f6f3", col = NA))
  # (a) input: genes x experiments heat strip computed from the toy X
  txt("Input: expression compendium", 0.2, Y(0.46), 7.5, face = "bold", just = "left")
  Xz <- t(scale(t(X))); Xz <- pmin(pmax(Xz, -2.5), 2.5); Xv <- (Xz + 2.5) / 5
  hx0 <- 0.3; hy <- 0.54; hw <- 1.0; hh <- 0.46
  cx <- rep(seq_len(N), each = G); cy <- rep(seq_len(G), times = N)
  rect(hx0 + (cx - 1) * hw / N, Y(hy) - cy * hh / G, hw / N + 0.002, hh / G + 0.002,
       fill = ramps$X[pmax(1, ceiling(as.vector(Xv) * 100))])
  txt("genes", hx0 - 0.04, Y(hy + hh / 2), 7, ink[["secondary"]], rot = 90, just = "right")
  txt("experiments", hx0 + hw / 2, Y(hy + hh + 0.1), 7, ink[["secondary"]])
  txt("toy: 14 x 60; real: thousands of genes,", 0.2, Y(hy + hh + 0.22), 7, ink[["secondary"]], just = "left")
  txt("hundreds+ experiments, many labs", 0.2, Y(hy + hh + 0.33), 7, ink[["secondary"]], just = "left")
  head_arrow(1.38, Y(0.77), 1.62, Y(0.77), ink[["secondary"]], 1.4)
  # (b) the four steps; published versus added bracket (thin grey = published, thick black = added)
  cg <- 0.09; cx0 <- 1.7; cy_ <- 0.95
  stepc <- list(list("Dependence", "MI (B-spline)", pal[["blue"]]),
                list("CLR", "z-score per gene", pal[["orange"]]),
                list("Softmax", "per-gene temp.", pal[["violet"]]),
                list("Walk", "t steps + residual", pal[["magenta"]]))
  tw <- function(l, face = "plain") convertWidth(grobWidth(textGrob(l, gp = gpar(fontfamily = FONT, fontsize = 7,
                                                                                fontface = face))), "native", TRUE)
  cws <- vapply(stepc, function(z) max(tw(z[[1]], "bold") + 0.28, tw(z[[2]]) + 0.1), 0)
  cxs <- cx0 + c(0, cumsum(cws[1:3] + cg))
  stopifnot(cxs[4] + cws[4] < 5.52)
  for (q in 1:4) {
    xq <- cxs[q]; cw <- cws[q]
    rect(xq, Y(cy_) - 0.2, cw, 0.4, fill = "white", col = ink[["grid"]], lwd = 0.8)
    rect(xq, Y(cy_) - 0.2, cw, 0.04, fill = stepc[[q]][[3]])
    chip(q, xq + 0.1, Y(cy_) + 0.09, 0.075)
    txt(stepc[[q]][[1]], xq + 0.2, Y(cy_) + 0.09, 7, face = "bold", just = "left")
    txt(stepc[[q]][[2]], xq + cw / 2, Y(cy_) - 0.06, 7, ink[["secondary"]])
    if (q < 4) head_arrow(xq + cw + 0.005, Y(cy_), xq + cw + cg - 0.005, Y(cy_), ink[["secondary"]], 1)
  }
  by_ <- Y(0.67)
  for (b in list(list(1, 2, ink[["secondary"]], 1, "published CLR (2007), modified", "plain"),
                 list(3, 4, ink[["primary"]], 2.4, "added here", "bold"))) {
    x0b <- cxs[b[[1]]]; x1b <- cxs[b[[2]]] + cws[b[[2]]]
    seg(x0b, by_, x1b, by_, b[[3]], b[[4]]); seg(x0b, by_, x0b, by_ - 0.04, b[[3]], b[[4]])
    seg(x1b, by_, x1b, by_ - 0.04, b[[3]], b[[4]])
    txt(b[[5]], (x0b + x1b) / 2, by_ + 0.09, 7, if (b[[6]] == "bold") ink[["primary"]] else ink[["secondary"]], b[[6]])
  }
  pill("no learned parameters  |  unsupervised  |  expression data only", (cx0 + cxs[4] + cws[4]) / 2, Y(1.28),
       h = 0.2, size = 7, fill = "white", col = ink[["primary"]], face = "bold", border = ink[["muted"]])
  head_arrow(cxs[4] + cws[4] + 0.06, Y(cy_), cxs[4] + cws[4] + 0.2, Y(cy_), ink[["secondary"]], 1.4)
  # (c) output: plain CLR S next to K_9 (each scaled to its own maximum), planted modules outlined
  ox <- c(5.62, 6.38); oc <- 0.033; oy <- 0.56; om <- G * oc
  txt("Output: gene x gene map", 5.6, Y(0.46), 7.5, face = "bold", just = "left")
  K9o <- Kt[["9"]]; diag(K9o) <- 0
  heat(ox[1], Y(oy), oc, S, ramps$S, vmax = max(S), gamma = 0.7)
  heat(ox[2], Y(oy), oc, K9o, ramps$K, vmax = max(K9o), gamma = 0.7)
  for (q in 1:2) for (blk in list(chain, modB))
    rect(ox[q] + (min(blk) - 1) * oc, Y(oy) - max(blk) * oc, length(blk) * oc, length(blk) * oc, col = "black", lwd = 0.9)
  txt("published CLR", ox[1] + om / 2, Y(oy + om + 0.09), 7, ink[["secondary"]])
  txt("this work", ox[2] + om / 2, Y(oy + om + 0.09), 7, ink[["secondary"]], "bold")
  txt("outlines: planted modules", 5.6, Y(oy + om + 0.21), 7, ink[["secondary"]], just = "left")
  txt("(e.g. regulons; not detected)", 5.6, Y(oy + om + 0.32), 7, ink[["secondary"]], just = "left")

  # ================= grid dividers ====================================================================
  r1 <- 1.5; r2 <- r1 + 2.8
  seg(0.1, Y(r2 - 0.07), W - 0.1, Y(r2 - 0.07), ink[["grid"]], 1)
  seg(XD, Y(r1), XD, Y(r2 - 0.07), ink[["grid"]], 1)
  seg(XD, Y(r2 - 0.07), XD, 0.32, ink[["grid"]], 1)

  # ============ ROW 1 =====================================================================
  yh <- Y(r1 + 0.1)
  # ---- panel 1: D -----------------------------------------------------------------------
  c1 <- 0.1; m1 <- G * c1; ytop1 <- Y(r1 + 0.8); x1 <- 0.72
  panel_head("1", "Dependence", 0.1, yh, "pairwise compatibility (parameter-free, symmetric)", pal[["blue"]],
             "MI (B-spline estimator here)")
  heat(x1, ytop1, c1, D, ramps$D, vmax = max(D))
  bands(x1, ytop1, c1, G, I, J, ext = 0.1)
  mark(x1, ytop1, c1, I, J)
  txt("gene i", x1 - 0.14, ytop1 - (I - 0.5) * c1, 7.5, face = "bold.italic", just = "right")
  txt("gene j", x1 + (J - 0.5) * c1, ytop1 + 0.1, 7.5, face = "bold.italic")
  txt("k", x1 + (K - 0.5) * c1, ytop1 + 0.1, 7.5, face = "bold.italic")
  tri(x1 + (K - 0.5) * c1, ytop1 + 0.03, 4)
  for (blk in list(chain, modB)) {                      # modules as brackets under the matrix
    xb0 <- x1 + (min(blk) - 1) * c1 + 0.01; xb1 <- x1 + max(blk) * c1 - 0.01; yb <- ytop1 - m1 - 0.15
    seg(xb0, yb, xb1, yb, ink[["primary"]], 1); seg(xb0, yb, xb0, yb + 0.04, ink[["primary"]], 1)
    seg(xb1, yb, xb1, yb + 0.04, ink[["primary"]], 1)
    txt("module", (xb0 + xb1) / 2, yb - 0.1, 7, ink[["secondary"]])
  }
  ptext(substitute(D[ij] == MI * "(" * i * ", " * j * ")" ~ "=" ~ v, list(v = sprintf("%.2f", D[I, J]))),
        x1 + m1 / 2, ytop1 - m1 - 0.42, 8)
  sw <- c("dependence measure:", "MI (B-spline) here;", "|r| and normal-", "scores correlation", "also tested")
  rect(2.3, ytop1 - 0.5, 0.03, 0.5, fill = pal[["blue"]])
  for (q in seq_along(sw)) txt(sw[q], 2.38, ytop1 - 0.06 - (q - 1) * 0.105, 7, ink[["secondary"]], just = "left")
  big_arrow(3.02, 3.52, ytop1 - m1 / 2 - 0.1, pal[["orange"]], h = 0.22, hl = 0.15)

  # ---- panel 2: CLR ------------------------------------------------------------------------
  panel_head("2", "Per-gene calibration (CLR)", XD + 0.1, yh, "per-query normalization of the logits", pal[["orange"]],
             c("CLR z-scores (Faith 2007);", "Stouffer combination here"))
  dx0 <- XD + 0.2; dx1 <- XD + 1.95; dh <- 0.3
  y_a <- Y(r1 + 1.22); y_b <- Y(r1 + 2.37)
  for (b in list(list(y_a, "i"), list(y_b, "j"))) {
    yb0 <- b[[1]]
    seg(dx0, yb0 + dh * 1.75 + 0.14, dx0 + 0.2, yb0 + dh * 1.75 + 0.14, "black", 1.3, lty = if (b[[2]] == "i") 1 else 2)
    ptext(substitute("row " * italic(g) * " of D", list(g = b[[2]])), dx0 + 0.25, yb0 + dh * 1.75 + 0.14, 7.5, just = "left")
  }
  di <- dist_panel(dx0, dx1, y_a, dh, zi, J, I)
  dj <- dist_panel(dx0, dx1, y_b, dh, zj, I, J)
  for (d_ in list(list(di, y_a, "i", zi), list(dj, y_b, "j", zj))) {
    dd <- d_[[1]]; yb0 <- d_[[2]]; g <- d_[[3]]; z <- d_[[4]]
    ptext(substitute(mu[g], list(g = as.name(g))), dd$X(dd$mu) - 0.03, yb0 + dh + 0.06, 7.5, just = "right")
    ptext(substitute(sigma[g], list(g = as.name(g))), dd$X(dd$mu + dd$s / 2), yb0 + dh * exp(-0.5) - 0.09, 7.5)
    ptext(substitute(zv ~ sigma[g], list(g = as.name(g), zv = sprintf("%.1f", z$z[if (g == "i") J else I]))),
          (dd$X(dd$mu) + dd$X(dd$obs)) / 2, dd$yb + 0.1, 7.5)
    ptext(expression(D[ij]), dd$X(dd$obs) - 0.04, yb0 + 0.17, 7.5, just = "right")
    txt("MI", dx0, yb0 - 0.12, 7, ink[["secondary"]], just = "left")
  }
  txt("MI < 0", di$X(LO) + (di$X(0) - di$X(LO)) / 2, y_a + dh * 1.02, 7, ink[["secondary"]])
  ptext(substitute(z[ij] == (D[ij] - mu[i]) / sigma[i] ~ "=" ~ v, list(v = sprintf("%.1f", zi$z[J]))),
        dx0, y_a - 0.29, 7.5, just = "left")
  ptext(substitute(z[ji] == (D[ij] - mu[j]) / sigma[j] ~ "=" ~ v, list(v = sprintf("%.1f", zj$z[I]))),
        dx0, y_b - 0.29, 7.5, just = "left")

  # S matrix (right column), the formula above it, two lines converging on the block arrow
  c2 <- 0.078; m2 <- G * c2; x2 <- XD + 2.3; ytop2m <- Y(r1 + 1.32)
  ptext(expression(S[ij] == frac(z[ij]^"+" + z[ji]^"+", sqrt(2))), x2 + m2 / 2, Y(r1 + 0.68), 8)
  ptext(substitute((a + b) / sqrt(2) == v, list(a = sprintf("%.1f", zi$z[J]), b = sprintf("%.1f", zj$z[I]),
                                               v = sprintf("%.1f", S[I, J]))), x2 + m2 / 2, Y(r1 + 0.95), 8)
  txt("z⁺ = max(z, 0)", x2 + m2 / 2, Y(r1 + 1.1), 7, ink[["secondary"]])
  cy <- ytop2m - m2 / 2
  seg(dx1 + 0.06, y_a + 0.0, XD + 2.17, cy + 0.04, ink[["secondary"]], 1.1)
  seg(dx1 + 0.06, y_b + 0.0, XD + 2.17, cy - 0.04, ink[["secondary"]], 1.1, lty = 2)
  big_arrow(XD + 2.17, x2 - 0.03, cy, pal[["orange"]], h = 0.2, hl = 0.1)
  heat(x2, ytop2m, c2, S, ramps$S, vmax = max(S))
  bands(x2, ytop2m, c2, G, I, J, ext = 0.06)
  mark(x2, ytop2m, c2, I, J, lwd = 1.5)
  txt("i", x2 - 0.09, ytop2m - (I - 0.5) * c2, 7.5, face = "bold.italic", just = "right")
  txt("j", x2 + (J - 0.5) * c2, ytop2m + 0.08, 7.5, face = "bold.italic")
  txt("S", x2 + m2 + 0.06, ytop2m - 0.05, 8, ink[["primary"]], "bold", just = "left")
  txt("white = clipped to 0", x2 + m2 / 2, ytop2m - m2 - 0.11, 7, ink[["secondary"]])
  txt("S is the input of step 3", x2 + m2 / 2, ytop2m - m2 - 0.23, 7, ORT, "bold")

  # ============ ROW 2 =====================================================================
  yh2 <- Y(r2 + 0.14); ytop3 <- Y(r2 + 0.82); c3 <- 0.068; gw <- G * c3
  # ---- panel 3: softmax -------------------------------------------------------------------
  panel_head("3", "Softmax attention", 0.1, yh2, "softmax over the keys, temperature per query", pal[["violet"]], "added here")
  xs3 <- 0.5; xa3 <- 2.15
  heat(xs3, ytop3, c3, S, ramps$S, vmax = max(S))
  bands(xs3, ytop3, c3, G, I, J, ext = 0.05)
  mark(xs3, ytop3, c3, I, J, lwd = 1.4)
  heat(xa3, ytop3, c3, A, ramps$A, vmax = max(A), gamma = 0.6)
  bands(xa3, ytop3, c3, G, I, J, ext = 0.05)
  mark(xa3, ytop3, c3, I, J, lwd = 1.4)
  for (x0_ in c(xs3, xa3)) {
    txt("i", x0_ - 0.09, ytop3 - (I - 0.5) * c3, 7.5, face = "bold.italic", just = "right")
    txt("j", x0_ + (J - 0.5) * c3, ytop3 + 0.08, 7.5, face = "bold.italic")
    txt("k", x0_ + (K - 0.5) * c3, ytop3 + 0.08, 7.5, face = "bold.italic")
  }
  txt("S, from step 2", xs3 + gw / 2, ytop3 + 0.21, 7.5, ORT, "bold")
  txt("A", xa3 + gw / 2, ytop3 + 0.21, 8, ink[["primary"]], "bold")
  ay <- ytop3 - 0.5 * gw
  big_arrow(xs3 + gw + 0.07, xa3 - 0.07, ay, pal[["violet"]], h = 0.2, hl = 0.1)
  txt("softmax", (xs3 + gw + xa3) / 2, ay + 0.17, 7, pal[["violet"]], "bold")
  # row i of S and of A as bar strips under the matrices, aligned to the columns
  by3 <- ytop3 - gw - 0.5; bh <- 0.28
  keep <- which(A[I, ] > 0)
  for (x0_ in c(xs3, xa3)) {
    rect(x0_ + (K - 1) * c3, by3, c3, bh, fill = ink[["muted"]], alpha = 0.3)                     # gene k slot
    rect(x0_ + (J - 1) * c3, by3, c3, bh, col = "black", lwd = 0.8, lty = 2)                      # gene j slot
    tri(x0_ + (K - 0.5) * c3, by3 - 0.08, 4)
  }
  sv <- S[I, ]; vm <- max(sv)
  for (g in seq_len(G)) if (sv[g] > 0)
    rect(xs3 + (g - 1) * c3 + 0.01, by3, c3 - 0.02, sv[g] / vm * bh,
         fill = if (g %in% keep) pal[["orange"]] else "#f6c9b2")
  seg(xs3, by3, xs3 + gw, by3, ink[["secondary"]], 0.8)
  for (g in seq_len(G)) if (A[I, g] > 0)
    rect(xa3 + (g - 1) * c3 + 0.01, by3, c3 - 0.02, A[I, g] / max(A[I, ]) * bh, fill = pal[["violet"]])
  seg(xa3, by3, xa3 + gw, by3, ink[["secondary"]], 0.8)
  seg(xa3 + (K - 1) * c3, by3, xa3 + K * c3, by3, ink[["primary"]], 2.4)       # A_ik = 0: hairline at the base
  txt("6th", xs3 + (K - 0.5) * c3 + 0.07, by3 - 0.08, 7, ink[["primary"]], "bold", just = "left")   # the cut bar is too short to see
  ptext(expression("row " * italic(i) * " of " * S), xs3 + gw / 2, by3 + bh + 0.1, 7.5)
  ptext(expression("row " * italic(i) * " of " * A * ", sums to 1"), xa3 + gw / 2, by3 + bh + 0.1, 7.5)
  ptext(substitute(S[ij] == v, list(v = sprintf("%.1f", S[I, J]))), xs3 + (J - 1) * c3 - 0.03,
        by3 + bh - 0.06, 7.5, just = "right")
  ptext(substitute(A[ij] == v, list(v = sprintf("%.2f", A[I, J]))), xa3 + (J - 1) * c3 - 0.03,
        by3 + bh - 0.06, 7.5, just = "right")
  txt("dark: top 5 kept (toy cap)", xs3 + gw / 2, by3 - 0.2, 7, ink[["secondary"]])
  txt("light: positive but cut", xs3 + gw / 2, by3 - 0.31, 7, ink[["secondary"]])
  ptext(substitute(paste(italic(k), " (triangle): ", S[ik] > 0), list()), xa3 + gw / 2, by3 - 0.2, 7)
  ptext(substitute(paste("6th, cut by cap: ", A[ik] == 0), list()), xa3 + gw / 2, by3 - 0.31, 7)
  ptext(substitute(paste(A[i] == softmax(S[i] / tau[i]), ";  ", tau[i], " set so that ", k[eff] == v),
                   list(v = KEFF)), (xs3 + xa3 + gw) / 2, by3 - 0.5, 7.5)
  ptext(expression("(toy " * k[eff] * " 3, cap 5; real method: " * k[eff] * " 10, top 50)"),
        (xs3 + xa3 + gw) / 2, by3 - 0.63, 7, ink[["secondary"]])

  # ---- panel 4: walk -----------------------------------------------------------------------
  c4 <- c3; gw4 <- G * c4; step4 <- gw4 + 0.13; x4 <- XD + 0.2
  big_arrow(xa3 + gw + 0.07, XD + 0.17, ytop3 - gw / 2, pal[["magenta"]], h = 0.2, hl = 0.1)
  panel_head("4", "Walk: depth", XD + 0.1, yh2, "repeated attention with residual", pal[["magenta"]], "added here")
  ptext(expression(P == 0.5 * I + 0.5 * A * ", repeated " * t * " times"), XD + 3.35, yh2 - 0.27, 7.5, just = "right")
  xs4 <- x4 + (0:2) * step4
  stopifnot(max(xs4) + gw4 < W - 0.05)
  offs <- lapply(Kt, function(m) { diag(m) <- 0; m }); kmax <- max(unlist(lapply(offs, function(m) m[I, ])))
  bbh <- 0.28; by4 <- ytop3 - gw4 - 0.2 - bbh
  for (q in seq_along(TS)) {
    M <- Kt[[q]]; off <- offs[[q]]
    heat(xs4[q], ytop3, c4, M, ramps$K, vmax = max(off), gamma = 0.7)
    bands(xs4[q], ytop3, c4, G, I, J, ext = 0.04, j_on = FALSE)
    rect(xs4[q] + (K - 1) * c4, ytop3 - I * c4, c4, c4, col = "black", lwd = 1.5)
    txt(sprintf("t = %d", TS[q]), xs4[q] + gw4 / 2, ytop3 + 0.17, 8.5, face = "bold")
    txt("k", xs4[q] + (K - 0.5) * c4, ytop3 + 0.05, 7, face = "bold.italic")
    if (TS[q] == TUSED)                      # framed: the depth the real method uses (the toy only displays it)
      rect(xs4[q] - 0.03, ytop3 - gw4 - 0.03, gw4 + 0.06, gw4 + 0.06 + 0.29, col = pal[["magenta"]], lwd = 2)
    rect(xs4[q] + (K - 1) * c4, by4, c4, bbh, fill = ink[["muted"]], alpha = 0.3)
    for (g in seq_len(G)) if (off[I, g] > 0)
      rect(xs4[q] + (g - 1) * c4 + 0.008, by4, c4 - 0.016, off[I, g] / kmax * bbh, fill = pal[["magenta"]])
    seg(xs4[q], by4, xs4[q] + gw4, by4, ink[["secondary"]], 0.8)
    if (off[I, K] == 0) seg(xs4[q] + (K - 1) * c4, by4, xs4[q] + K * c4, by4, ink[["primary"]], 2.4)
    tri(xs4[q] + (K - 0.5) * c4, by4 - 0.08, 4)
    ptext(substitute(K[ik] == v, list(v = sprintf("%.3f", M[I, K]))), xs4[q] + gw4 / 2, by4 - 0.2, 7.5)
    if (q < 3) head_arrow(xs4[q] + gw4 + 0.012, ytop3 - gw4 / 2, xs4[q + 1] - 0.012, ytop3 - gw4 / 2,
                          pal[["magenta"]], 1.1)
  }
  txt("i", xs4[1] - 0.06, ytop3 - (I - 0.5) * c4, 7.5, face = "bold.italic", just = "right")
  ptext(expression("row " * italic(i) * " of " * K[t] == (P^t + (P^t)^T) / 2 * ";  triangle: gene " * italic(k)),
        xs4[1], by4 + bbh + 0.1, 7, ink[["secondary"]], just = "left")
  # what the walk buys: co-membership AP, plain CLR versus K_t (illustrative, favourable by construction)
  ap0 <- by4 - 0.66; apw <- gw4 * 0.72
  txt("toy: average precision (AP) of ranking same-module pairs;", xs4[1], by4 - 0.36, 7, ink[["secondary"]],
      "bold", just = "left")
  txt("favourable by construction. Real data and losses: Figs. 2-6", xs4[1], by4 - 0.47, 7, ink[["secondary"]],
      "bold", just = "left")
  for (q in seq_along(TS)) {
    rect(xs4[q], ap0 - 0.04, AP_K[[q]] * apw, 0.08, fill = pal[["magenta"]])
    seg(xs4[q] + AP_BASE * apw, ap0 - 0.05, xs4[q] + AP_BASE * apw, ap0 + 0.05, ink[["muted"]], 1)
    seg(xs4[q] + AP_CLR * apw, ap0 - 0.08, xs4[q] + AP_CLR * apw, ap0 + 0.08, ORT, 2.4)
    txt(sprintf("%.2f", AP_K[[q]]), xs4[q] + apw + 0.05, ap0, 7, just = "left")
  }
  txt(sprintf("rand. %.2f", AP_BASE), xs4[1], ap0 - 0.14, 7, ink[["secondary"]], just = "left")
  txt(sprintf("CLR %.2f", AP_CLR), xs4[1] + AP_CLR * apw - 0.03, ap0 - 0.14, 7, ORT, "bold", just = "left")
  txt("framed: depth used by the method (chosen without labels on", xs4[1], ap0 - 0.29, 7, ink[["secondary"]], just = "left")
  txt("held-out experiments); the toy only displays it", xs4[1], ap0 - 0.39, 7, ink[["secondary"]], just = "left")

  # ============ footer ================================================================================
  fy <- 0.3
  seg(0.1, fy, W - 0.1, fy, ink[["grid"]], 1)
  txt("Toy example; analogy with transformer attention, not equivalence; no real-data result shown.",
      0.1, fy - 0.14, 7, ink[["primary"]], just = "left")
  popViewport()
}

dev_png <- function() png(file.path(out_dir, "fig1_method.png"), width = W * 300, height = H * 300,
                          res = 300, type = "cairo", bg = "white")
dev_png(); draw_fig(); invisible(dev.off())
cairo_pdf(file.path(out_dir, "fig1_method.pdf"), width = W, height = H, bg = "white")
draw_fig(); invisible(dev.off())
cat("wrote", file.path(out_dir, "fig1_method.{png,pdf}"), "; height", H, "in\n")
