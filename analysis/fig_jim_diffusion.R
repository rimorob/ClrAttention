# Attention diffusing over a gene network (multi-step lazy walk), for the PI's private
# readout to J. Collins ONLY; not for submission. Drawn in code with grid; layout and
# affinities are a seeded SYNTHETIC toy (no data), but every mass value drawn is computed
# with the method's own rules: per-gene softmax over the <= 50 strongest candidates with the
# temperature set so that exp(entropy) = 10 effective partners, lazy walk P = 0.5 I + 0.5 A,
# read-out = symmetrised mass (P^t + P^t')/2 for the query row. No learned parameters.
# Panel layout was guided by PaperBanana drafts (docs/whitepapers/figures/pb_runs/, untracked),
# none of whose content is used.
# Palette: analysis/plot_theme.R (magenta = depth/walk, grey = neutral).
#
# Output: docs/whitepapers/figures/jim_diffusion.png
# Run:    nice -n 19 Rscript analysis/fig_jim_diffusion.R

suppressPackageStartupMessages(library(grid))
source("analysis/plot_theme.R")
out_dir <- "docs/whitepapers/figures"
FONT <- "Liberation Sans"

# ---- method (matches the package rules) -------------------------------------
eff_n <- 10; max_cand <- 50
attention_matrix <- function(S) {          # S: symmetric relatedness score matrix (CLR-like, >= 0)
  n <- nrow(S); A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    s <- S[i, ]; s[i] <- -Inf
    cand <- order(s, decreasing = TRUE)[seq_len(min(max_cand, n - 1))]
    z <- s[cand]
    f <- function(lt) { w <- exp((z - max(z)) / exp(lt)); p <- w / sum(w); p <- p[p > 0]; exp(-sum(p * log(p))) - eff_n }
    lt <- uniroot(f, c(-5, 3))$root                         # per-gene temperature
    w <- exp((z - max(z)) / exp(lt)); A[i, cand] <- w / sum(w)
  }
  A
}
lazy_walk <- function(A, i, steps) {       # symmetrised mass of query row after each step count
  n <- nrow(A); P <- 0.5 * diag(n) + 0.5 * A; Pt <- diag(n); out <- list(); k <- 0
  for (t in 0:max(steps)) {
    if (t > 0) Pt <- Pt %*% P
    if (t %in% steps) out[[as.character(t)]] <- (Pt[i, ] + Pt[, i]) / 2
  }
  out
}

# ---- toy network 1: a loose regulon among neighbouring modules ----------------
set.seed(11)
mk_module <- function(n, cx, cy, sd, dmin = 0.22) {        # Gaussian cloud with a minimum spacing
  P <- NULL
  while (is.null(P) || nrow(P) < n) {
    q <- c(rnorm(1, cx, sd), rnorm(1, cy, sd))
    if (is.null(P) || min(sqrt((P[, 1] - q[1])^2 + (P[, 2] - q[2])^2)) > dmin) P <- rbind(P, q)
  }
  P
}
mods <- list(reg = mk_module(18, 0, 0, 0.62, 0.30), m1 = mk_module(11, 1.95, 0.55, 0.34),
             m2 = mk_module(10, -1.85, 0.7, 0.32), m3 = mk_module(12, 0.2, -1.95, 0.36),
             m4 = mk_module(10, -0.5, 2.0, 0.30), m5 = mk_module(10, 1.75, -1.55, 0.30))
X <- do.call(rbind, mods); mod_id <- rep(names(mods), sapply(mods, nrow))
X <- X[, 1:2]
# query gene i: the regulon member closest to the regulon centre
reg_idx <- which(mod_id == "reg")
i0 <- reg_idx[which.min(rowSums(X[reg_idx, ]^2))]
D2 <- as.matrix(dist(X))^2
S1 <- exp(-D2 / (2 * 0.75^2)) * 6          # CLR-like score: decays with distance
set.seed(12); N <- matrix(rnorm(nrow(X)^2, 0, 0.35), nrow(X)); S1 <- pmax(S1 + (N + t(N)) / 2, 0)
A1 <- attention_matrix(S1)
mass1 <- lazy_walk(A1, i0, c(0, 1, 3, 9))

cat("network 1: regulon share of symmetrised mass at t = 0,1,3,9:",
    sapply(mass1, function(m) round(sum(m[reg_idx]) / sum(m), 2)), "\n")

# ---- toy network 2: a tight clique (complex) inside a wider program ------------
set.seed(21)
nP <- 40; nC <- 6
ang <- seq(0, 2 * pi, length.out = nC + 1)[-(nC + 1)] + pi / 2
XC <- cbind(-0.55 + 0.42 * cos(ang), -0.15 + 0.42 * sin(ang))
XP <- cbind(rnorm(nP, 0.15, 0.95), rnorm(nP, 0.1, 0.75))
keep <- apply(XP, 1, function(p) min(sqrt((p[1] - XC[, 1])^2 + (p[2] - XC[, 2])^2)) > 0.38)
XP <- XP[keep, ]; nP <- nrow(XP)
X2 <- rbind(XC, XP); grp2 <- c(rep("clique", nC), rep("program", nP))
Dd <- as.matrix(dist(X2))
S2 <- 3 * exp(-Dd^2 / (2 * 0.8^2))                    # program: score decays with distance
S2[1:nC, 1:nC] <- 6                                    # tight complex: all pairs strongly related
set.seed(22); N <- matrix(rnorm(nrow(X2)^2, 0, 0.3), nrow(X2)); S2 <- pmax(S2 + (N + t(N)) / 2, 0)
A2 <- attention_matrix(S2)
q2 <- 1
mass2 <- lazy_walk(A2, q2, c(0, 9))
cat("network 2: clique share of mass at t = 0, 9:",
    sapply(mass2, function(m) round(sum(m[1:nC]) / sum(m), 2)), "\n")

# ---- drawing helpers ----------------------------------------------------------
mag <- pal[["magenta"]]; mag_deep <- "#b23a72"
ramp <- colorRampPalette(c("#ffffff", "#f7cddd", mag, mag_deep))
mass_col <- function(m, mref) ramp(256)[1 + round(255 * pmin(1, (m / mref)^0.45))]   # compressive scale: thin tails stay visible
txt <- function(label, x, y, size = 14, col = ink[["primary"]], face = "plain", just = "centre")
  grid.text(label, x, y, default.units = "native", just = just,
            gp = gpar(fontfamily = FONT, fontsize = size, col = col, fontface = face))
edges_top <- function(A, k = 3, thr = 0.07) {   # strongest k partners per gene, plus any link with weight > thr
  e <- NULL
  for (i in seq_len(nrow(A))) { j <- order(A[i, ], decreasing = TRUE)[1:k]; e <- rbind(e, cbind(pmin(i, j), pmax(i, j))) }
  w <- which(A > thr, arr.ind = TRUE); e <- rbind(e, cbind(pmin(w[, 1], w[, 2]), pmax(w[, 1], w[, 2])))
  unique(e)
}
hull_shape <- function(P, pad = 0.28) {    # padded convex hull, smoothed
  h <- P[chull(P), , drop = FALSE]; ctr <- colMeans(h)
  d <- sweep(h, 2, ctr); h <- h + pad * d / sqrt(rowSums(d^2))
  h
}
glow <- function(P, w, xr, yr, bw = 0.30, n = 180, pad = 0.9) {   # field padded so no box edge shows
  xe <- xr + c(-pad, pad); ye <- yr + c(-pad, pad)
  gx <- seq(xe[1], xe[2], length.out = n); gy <- seq(ye[2], ye[1], length.out = n)
  F <- matrix(0, n, n)
  for (k in which(w > 1e-4 * max(w))) {
    kx <- exp(-(gx - P[k, 1])^2 / (2 * bw^2)); ky <- exp(-(gy - P[k, 2])^2 / (2 * bw^2))
    F <- F + w[k] * outer(ky, kx)
  }
  g <- (F / max(F))^0.7
  matrix(paste0(mag, sprintf("%02X", round(255 * 0.62 * g))), n, n)
}
# Draw one network in a viewport of the requested size (inches) centred at (cx, cy) inches.
draw_net <- function(cx, cy, wd, ht, P, edges, mass, mref, xr, yr, query = NULL, hull = NULL,
                     faded = NULL, hide = NULL, clique = NULL, leaks = NULL, rad = 0.085) {
  pushViewport(viewport(x = unit(cx, "in"), y = unit(cy, "in"), width = unit(wd, "in"),
                        height = unit(ht, "in"), xscale = xr, yscale = yr, clip = "off"))
  vis <- if (is.null(faded)) rep(1, nrow(P)) else ifelse(faded, 0.28, 1)
  if (!is.null(hide)) vis[hide] <- 0
  padx <- 0.9 / diff(xr); pady <- 0.9 / diff(yr)
  grid.raster(glow(P, (mass / mref)^0.45 * vis, xr, yr), x = 0.5, y = 0.5, width = 1 + 2 * padx,
              height = 1 + 2 * pady, interpolate = TRUE)
  if (!is.null(hull)) grid.xspline(hull[, 1], hull[, 2], default.units = "native", shape = 1, open = FALSE,
                                   gp = gpar(col = ink[["muted"]], lty = "dashed", lwd = 1.1, fill = NA))
  grid.segments(P[edges[, 1], 1], P[edges[, 1], 2], P[edges[, 2], 1], P[edges[, 2], 2],
                default.units = "native",
                gp = gpar(col = ink[["muted"]], lwd = 0.7, alpha = 0.6 * pmin(vis[edges[, 1]], vis[edges[, 2]])))
  if (!is.null(clique)) {
    cp <- t(combn(clique, 2))
    grid.segments(P[cp[, 1], 1], P[cp[, 1], 2], P[cp[, 2], 1], P[cp[, 2], 2], default.units = "native",
                  gp = gpar(col = ink[["secondary"]], lwd = 1.1, alpha = 0.85))
  }
  fills <- mass_col(mass, mref)
  ord <- seq_len(nrow(P))
  grid.circle(P[ord, 1], P[ord, 2], rad, default.units = "native",
              gp = gpar(fill = fills[ord], col = ink[["secondary"]], lwd = 0.9, alpha = vis))
  if (!is.null(leaks))
    for (r in seq_len(nrow(leaks)))
      grid.curve(leaks[r, 1], leaks[r, 2], leaks[r, 3], leaks[r, 4], default.units = "native",
                 curvature = 0.3, square = FALSE, ncp = 8, shape = 0.6,
                 arrow = arrow(length = unit(0.1, "in"), type = "closed"),
                 gp = gpar(col = mag_deep, fill = mag_deep, lwd = 2))
  if (!is.null(query)) {
    grid.circle(P[query, 1], P[query, 2], rad * 1.55, default.units = "native",
                gp = gpar(fill = fills[query], col = ink[["primary"]], lwd = 2.4))
    txt("i", P[query, 1] + rad * 3.0, P[query, 2] + rad * 3.2, 15, ink[["primary"]], "bold.italic")
  }
  popViewport()
}
panel_arrow <- function(x, y)
  grid.polygon(unit(x + c(-0.07, 0.07, -0.07), "in"), unit(y + c(-0.14, 0, 0.14), "in"),
               gp = gpar(fill = ink[["muted"]], col = NA))
rbox <- function(x0, y0, x1, y1)
  grid.roundrect(unit((x0 + x1) / 2, "in"), unit((y0 + y1) / 2, "in"), unit(x1 - x0, "in"), unit(y1 - y0, "in"),
                 r = unit(0.12, "in"), gp = gpar(fill = "white", col = "#d9d8d3", lwd = 1.2))

# ---- compose (16 x 9 in canvas, 2400 x 1350 px) -------------------------------
W <- 16; H <- 9
png(file.path(out_dir, "jim_diffusion.png"), width = 2400, height = 1350, res = 150, bg = ink[["surface"]])
grid.newpage(); pushViewport(viewport(xscale = c(0, W), yscale = c(0, H), width = 1, height = 1))
options_unit <- function(x, y) NULL
# native units here are inches because the root xscale is 0..W (16 in wide)
xr1 <- range(X[, 1]) + c(-0.35, 0.35); yr1 <- range(X[, 2]) + c(-0.35, 0.35)
asp1 <- diff(xr1) / diff(yr1)
pw <- 3.6; ph <- pw / asp1
steps <- c("0", "1", "3", "9"); lab <- paste("step", steps)
e1 <- edges_top(A1, 3)
hull1 <- hull_shape(X[reg_idx, ], 0.30)
pitch <- 3.95; x0 <- 0.45 + pw / 2; ytop <- 6.2
for (k in seq_along(steps)) {
  cx <- x0 + (k - 1) * pitch
  m <- mass1[[steps[k]]]
  draw_net(cx, ytop, pw, ph, X, e1, m, max(m), xr1, yr1, query = i0, hull = hull1)
  txt(lab[k], cx, ytop + ph / 2 + 0.30, 17, ink[["primary"]], "bold")
  if (k < 4) panel_arrow(cx + pitch / 2, ytop)
}
# ---- bottom row: the trade-off, both after 9 steps ------------------------------
by0 <- 0.50; by1 <- 3.95
rbox(0.4, by0, 7.85, by1); rbox(8.15, by0, 15.6, by1)
txt("regulon (loose)", 0.65, by1 - 0.32, 16, ink[["primary"]], "bold", just = "left")
txt("complex (tight clique)", 8.40, by1 - 0.32, 16, ink[["primary"]], "bold", just = "left")
txt("step 9", 7.60, by1 - 0.32, 12.5, ink[["secondary"]], just = "right")
txt("step 9", 15.35, by1 - 0.32, 12.5, ink[["secondary"]], just = "right")
# left: the regulon, zoomed; neighbours faded
m9 <- mass1[["9"]]
xrl <- range(X[reg_idx, 1]) + c(-0.5, 0.5); yrl <- range(X[reg_idx, 2]) + c(-0.45, 0.45)
inl <- which(X[, 1] > xrl[1] & X[, 1] < xrl[2] & X[, 2] > yrl[1] & X[, 2] < yrl[2])
hl <- 2.6; wl <- hl * diff(xrl) / diff(yrl)
edl <- e1[e1[, 1] %in% inl & e1[, 2] %in% inl, , drop = FALSE]
draw_net(4.12, 1.95, wl, hl, X, edl, m9, max(m9[reg_idx]), xrl, yrl, query = i0, hull = hull1,
         faded = mod_id != "reg", hide = setdiff(seq_len(nrow(X)), inl), rad = 0.07)
# right: the clique inside the wider program
m92 <- mass2[["9"]]
xr2 <- range(X2[, 1]) + c(-0.25, 0.25); yr2 <- range(X2[, 2]) + c(-0.25, 0.25)
h2 <- 2.6; w2 <- h2 * diff(xr2) / diff(yr2)
e2 <- edges_top(A2, 3); e2 <- e2[!(e2[, 1] <= nC & e2[, 2] <= nC), , drop = FALSE]
hc <- hull_shape(X2[1:nC, ], 0.22)
cc <- colMeans(X2[1:nC, ])
far <- nC + which(sqrt(rowSums(sweep(X2[-(1:nC), ], 2, cc)^2)) > 1.0)     # program genes well outside the clique
far <- far[order(m92[far], decreasing = TRUE)]
tg <- far[1]                                                        # three highest-mass genes in distinct directions
for (k in far[-1]) {
  a <- atan2(X2[k, 2] - cc[2], X2[k, 1] - cc[1])
  if (all(abs(sin((a - sapply(tg, function(q) atan2(X2[q, 2] - cc[2], X2[q, 1] - cc[1]))) / 2)) > 0.45)) tg <- c(tg, k)
  if (length(tg) == 3) break
}
lk <- t(sapply(tg, function(k) { d <- X2[k, ] - cc; u <- d / sqrt(sum(d^2))
  c(cc + 0.40 * u, X2[k, ] - 0.14 * u) }))
draw_net(11.9, 1.95, w2, h2, X2, e2, m92, max(m92[1:nC]), xr2, yr2, query = q2, hull = hc,
         clique = 1:nC, leaks = lk, rad = 0.07)
txt("schematic", 0.40, 0.22, 11, ink[["muted"]], "italic", just = "left")
dev.off()
