# Library for the exploratory tcaCLR depth / width / clipping grid (docs/plans/tcaclr_depth_grid_plan.md).
# EXPLORATORY; nothing here is registered and the registered scorer (analysis/tcaclr_confirm/*.R) is not modified.
# Source from the repo root:  source("analysis/tcaclr_depth/grid_lib.R")
#
# Cell = (method aCLR|tcaCLR, keff, cap, clipS, clipO, t).  An OPERATOR = (keff, cap, clipS); it carries all depths and both clipO.
#   clipS  TRUE : pairwise score S = (max(0, z_ij) + max(0, z_ji)) / sqrt(2)   (= the registered CLR score; operator built exactly as
#                 ClrAttention$build_operator(softmax_keff, softmax_cap) does: clr:::.softmax_rows, row normalisation, self-loop)
#          FALSE: S = (z_ij + z_ji) / sqrt(2) from the same calibrated rows without clipping; candidates are still the top-cap partners with S > 0
#   clipO  TRUE : tcaCLR output Z clipped at 0 before the Stouffer combination (registered: sk_a3r)
#          FALSE: Z kept negative (NA = unreachable pairs still set to 0, as in the registered read-out)
# The clipped path calls only registered code (clr:::.softmax_rows, the finish step copied line by line, the lazy walk of make_scores.R,
# sk_a3r arithmetic) and is checked for bit-identity in check_reference.R.

suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R")   # ap_fast, the registered tie convention

GL <- list(
  # Revision 2 of the plan: 17 depths (t = 9 on the grid), 8 widths, cap = 4 x keff (at most 320); the reference cell keeps cap 50
  depths = c(1L, 2L, 3L, 4L, 6L, 8L, 9L, 12L, 13L, 16L, 18L, 24L, 32L, 48L, 64L, 96L, 128L),
  keffs  = c(7, 10, 14, 20, 28, 40, 56, 80),
  cap_mult = 4L, cap_max = 320L, alpha = 0.5,
  eps = c(0.01, 0.1),                           # mixture floor of C3: p_mix = (1 - eps) p + eps / (G - 1); eps is a parameter, two values reported
  bic_k = 3L,                                   # hyperparameters counted in the BIC reading of C3 (keff, t, clip); never combined with C3
  ref = list(keff = 10, cap = 50L, t = 9L),     # the registered cell: keff 10, cap 50, t 9, clipped both
  top_frac = 0.01)                              # C1/C2 target: top 1% of the other half's score

## ---- fitting ------------------------------------------------------------------------------------------------------------------
gl_fit <- function(X, threads = 8L) {            # as make_scores.R / lab_splits.R build_fit
  ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none", threads = threads)$calibrate(method = "normal", combine = "stouffer")
}
# S = (z_ij + z_ji)/sqrt(2), z = row z-scores of the MI matrix with zeroed diagonal (mean and sd over all G entries, sd with G - 1), no clipping.
# Same formula as src/core/clr_core.cpp clr_calibrate with the clipping line removed (and the same multiplication by 1/sqrt(2)).
gl_unclipped_clr <- function(M) {
  G <- nrow(M); diag(M) <- 0
  mu <- rowSums(M) / G; M <- M - mu
  sdv <- sqrt(rowSums(M^2) / (G - 1))
  z <- M / sdv; z[!is.finite(z)] <- 0; rm(M)
  S <- (z + t(z)) * (1 / sqrt(2)); diag(S) <- 0
  S
}
gl_pair_score <- function(f, clipS) if (clipS) f$clr_scores else gl_unclipped_clr(f$mi)

## ---- operator -----------------------------------------------------------------------------------------------------------------
# finish step of ClrAttention$.finish_operator (row-normalise, rows without mass attend to themselves), line by line
gl_finish_A <- function(A) {
  rs <- rowSums(A); nz <- rs > 0
  A[nz, ] <- A[nz, ] / rs[nz]
  if (any(!nz)) A[cbind(which(!nz), which(!nz))] <- 1
  A
}
gl_operator_A <- function(S, keff, cap) gl_finish_A(clr:::.softmax_rows(S, keff = keff, cap = as.integer(cap)))
# lazy operator P = (1 - alpha) I + alpha A, sparse, as make_scores.R
gl_make_P <- function(S, keff, cap, alpha = GL$alpha) {
  A <- gl_operator_A(S, keff, cap); G <- nrow(A)
  P <- Matrix::Diagonal(G, 1 - alpha) + alpha * Matrix::Matrix(A, sparse = TRUE); rm(A); P
}

# stationary distribution of the lazy walk with teleport 1e-3, as rz_stationary in analysis/robust_z/robust_z_lib.R (copied, same iteration)
gl_stationary <- function(P, eps = 1e-3) {
  G <- nrow(P); pi <- rep(1 / G, G)
  for (it in 1:20000) { q <- (1 - eps) * as.numeric(pi %*% P) + eps / G; if (sum(abs(q - pi)) < 1e-13) break; pi <- q }
  q
}

## ---- powers of P by the fewest dense products ---------------------------------------------------------------------------------
# Every depth t > 1 is computed from one earlier power: t even from (t/2)^2 (dense x dense), t odd from (t-1) %*% P (dense x sparse).
# The depths of GL$depths are closed under halving (4 = 2^2, 6 = 3^2, 12 = 6^2, ..., 256 = 128^2): 12 dense products plus the sparse ones
# for t = 2 and 3; a registered-style sequential chain would need 255. Extra depths (e.g. 9) pull in the powers they need (8, 4, 2).
# fn(t, Pt) is called for each requested depth, in increasing order, with the dense base matrix (read-only); a matrix is dropped as soon as
# nothing later needs it.
gl_dep <- function(t) if (t == 1L) NA_integer_ else if (t %% 2L == 0L) t %/% 2L else t - 1L
gl_need <- function(depths) {
  need <- sort(unique(as.integer(depths)))
  repeat { d <- stats::na.omit(vapply(need, gl_dep, 0L)); add <- setdiff(d, need); if (!length(add)) break; need <- sort(c(need, add)) }
  need
}
gl_n_dense <- function(depths) { need <- gl_need(depths); sum(need > 3L & need %% 2L == 0L) }   # dense x dense products
gl_walk <- function(P, depths, fn) {
  depths <- sort(unique(as.integer(depths))); need <- gl_need(depths)
  parent <- vapply(need, gl_dep, 0L); lastuse <- stats::setNames(rep(0L, length(need)), need)
  for (k in seq_along(need)) if (!is.na(parent[k])) lastuse[as.character(parent[k])] <- max(lastuse[as.character(parent[k])], need[k])
  mats <- list()
  for (t in need) {
    key <- as.character(t)
    m <- if (t == 1L) as.matrix(P) else {
      md <- mats[[as.character(gl_dep(t))]]
      if (t %% 2L == 0L) md %*% md else as.matrix(md %*% P) }
    mats[[key]] <- m; rm(m)
    if (t %in% depths) fn(t, mats[[key]])
    d <- gl_dep(t)
    if (!is.na(d) && lastuse[[as.character(d)]] <= t) mats[[as.character(d)]] <- NULL
    if (lastuse[[key]] == 0L) mats[[key]] <- NULL
  }
  invisible(NULL)
}
# the registered chain (make_scores.R): t - 1 sparse products one after the other (used by the reference check)
gl_power_seq <- function(P, t) { Pt <- as.matrix(P); for (s in seq_len(t - 1L)) Pt <- as.matrix(Pt %*% P); Pt }

## ---- score read-outs ----------------------------------------------------------------------------------------------------------
# upper-triangle index (column-major, = upper.tri order = tc_pack) and its transpose; ri < cj are the pair genes
gl_ut <- function(G) {
  G <- as.integer(G); cj <- rep.int(seq_len(G)[-1L], times = seq_len(G - 1L)); ri <- sequence(seq_len(G - 1L))
  list(G = G, ut = ri + (cj - 1L) * G, utT = cj + (ri - 1L) * G, ri = ri, cj = cj)
}
# Row z of log P^t over reachable off-diagonal entries; unreachable (P = 0) and the diagonal are NA -> 0. Negative z are kept here
# (clipO = FALSE); clip with pmax(z, 0) for clipO = TRUE. Same operations in the same order as stack_lib.R sk_a3r.
gl_zmat <- function(Pt) {
  L <- log(Pt); L[!is.finite(L)] <- NA_real_; diag(L) <- NA_real_
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  L <- L - mu
  sdv <- sqrt(rowSums(L^2, na.rm = TRUE) / (n - 1))
  L <- L / sdv; L[is.na(L)] <- 0
  L
}
# upper-triangle score vectors of one P^t: aCLR = (P^t + t(P^t))/2 ; tcaCLR = (Z + t(Z))/sqrt(2) with Z clipped (clipO) or not
gl_variants <- function(Pt, U) {
  a <- (Pt[U$ut] + Pt[U$utT]) / 2
  Z <- gl_zmat(Pt); z1 <- Z[U$ut]; z2 <- Z[U$utT]; rm(Z)
  list(aCLR = a, tca_clip = (pmax(z1, 0) + pmax(z2, 0)) / sqrt(2), tca_unclip = (z1 + z2) / sqrt(2))
}
gl_variant_names <- c(aCLR = "aCLR", tca_clip = "tcaCLR", tca_unclip = "tcaCLR")
gl_variant_clipO <- c(aCLR = NA, tca_clip = TRUE, tca_unclip = FALSE)

## ---- criteria (plan Revision 2) ----------------------------------------------------------------------------------------------
# AP of score v (sorted ascending copy sv) for the pair set Tidx, with the registered tie convention (ties at the end of their block),
# computed from the positives only: identical arithmetic to ap_fast(v, label) for label = 1 on Tidx (tested bit-identical).
gl_ap_T <- function(v, sv, Tidx) {
  N <- length(v); m <- length(Tidx)
  vt <- sort.int(v[Tidx], decreasing = TRUE, method = "radix")
  chg <- c(vt[-1L] != vt[-m], TRUE); ends <- which(chg)
  tp <- as.numeric(ends); vd <- vt[ends]
  k <- N - findInterval(vd, sv, left.open = TRUE)
  sum(diff(c(0, tp)) * tp / k) / m
}
gl_sorted <- function(v) sort.int(v, method = "radix")
# top 1% pair set of a score (the convention of lab_splits.R / sk_top1): v >= quantile(v, 0.99), and v > 0 when positive_only
gl_top <- function(v, sv, frac = GL$top_frac, positive_only = TRUE) {
  q <- stats::quantile(sv, 1 - frac, names = FALSE)
  which(v >= q & (if (positive_only) v > 0 else TRUE))
}
# C1 (DIAGNOSTIC ONLY, plan Revision 2): AP of half 1's score for the top 1% of half 2's CLR (similarity to CLR; favours shallow settings).
# C2 (null-adjusted): AP of half 1's score for the top 1% of half 2's score by the same method, minus the same statistic for the
#   degree-only surrogate of the same operator (stationary lift: P^inf[i, j] = pi_j, read out by the same aCLR / tcaCLR formulas,
#   pi = stationary distribution of that half's lazy walk). gl_stage1_split returns C2 (raw), C2_null and C2_adj = C2 - C2_null.
# C3 (blocked held-out likelihood, no penalty): p[i, j] = P^t[i, j] / (sum over the reachable off-diagonal entries of row i) for both directed
#   instances (i -> j) and (j -> i) of every held-out target pair (top 1% of the other half's CLR); mixture floor
#   p_mix = (1 - eps) p + eps / (G - 1) (so an unreachable instance, p = 0, scores log(eps / (G - 1))); row-balanced weights: the mean log p_mix
#   within each row that has instances, then the mean over those rows. "Blocked": the halves are unions of laboratories (M3D) or knockdown
#   clusters (K562) and each direction is one fold. Reported at eps 0.01 (C3) and 0.1 (C3_e10); c3_unreach = share of unreachable instances.
# C3_bic (second reading only, never combined with C3): n_eff * C3 - 0.5 * k * log(n_eff), n_eff = effective number of independent experiments of the
#   held-out half ((tr C)^2 / tr(C^2) of the experiment correlation matrix), k = GL$bic_k hyperparameters. (k_reach, the mean reachable entries
#   per row of P^t, is a mixing diagnostic only.)
gl_c3 <- function(Pt, ti, tj) {
  G <- nrow(Pt); d <- diag(Pt)
  rs <- rowSums(Pt) - d
  k <- mean(rowSums(Pt > 0) - (d > 0))
  p <- c(Pt[ti + (tj - 1L) * G] / rs[ti], Pt[tj + (ti - 1L) * G] / rs[tj]); p[!is.finite(p)] <- 0
  rw <- c(ti, tj); cnt <- tabulate(rw, G); use <- which(cnt > 0)
  bal <- function(eps) { lp <- log((1 - eps) * p + eps / (G - 1)); mean((rowsum(lp, rw, reorder = FALSE)[, 1][as.character(use)] / cnt[use])) }
  c(C3 = bal(GL$eps[1]), C3_e10 = bal(GL$eps[2]), c3_unreach = mean(!(p > 0)), k_reach = k)
}
# effective sample sizes of an expression matrix X (genes x experiments): N_eff = (tr C)^2 / tr(C^2), C the experiment correlation matrix
# of gene-centred data (columns standardised over genes); G_eff the same for the gene correlation matrix (rows standardised over experiments). tr(C^2)
# is computed from the small experiments x experiments Gram matrix, so no G x G matrix is formed. G_eff is on the expression data, not on CLR rows.
gl_neff <- function(X) {
  n <- ncol(X); G <- nrow(X)
  Xc <- X - rowMeans(X)                                                # gene means removed first (otherwise every experiment correlates with the mean profile)
  Zc <- scale(Xc); Zc[!is.finite(Zc)] <- 0                             # columns standardised over genes
  C <- crossprod(Zc) / (G - 1); n_eff <- n^2 / sum(C^2)
  Zr <- t(scale(t(X))); Zr[!is.finite(Zr)] <- 0                        # rows standardised over experiments
  Cr <- crossprod(Zr) / (n - 1)                                        # n x n; ||Zr Zr'||_F^2 = ||Zr' Zr||_F^2
  g_eff <- G^2 / sum(Cr^2)                                             # = G^2 / ||gene correlation matrix||_F^2
  c(N = n, N_eff = n_eff, G = G, G_eff = g_eff)
}
gl_ap <- function(score, label) ap_fast(score, label)   # pooled AUPR on a truth: the registered convention

## ---- stage 1 for one split ---------------------------------------------------------------------------------------------------
# All operators of `ops` (list of list(keff, cap, clipS, depths)) for the two half fits f1, f2. Direction 1: half 1's score against half 2's
# targets; direction 2: the other way. n_eff = c(N_eff of half 1, N_eff of half 2) for the BIC reading. Returns one data.frame.
# Memory: the second half's score vectors of one operator (3 variants x depths) are held while the first half is walked.
gl_stage1_split <- function(f1, f2, ops, n_eff = c(NA_real_, NA_real_), say = function(...) invisible(NULL)) {
  G <- nrow(f1$clr_scores); U <- gl_ut(G)
  clr1 <- f1$clr_scores[U$ut]; clr2 <- f2$clr_scores[U$ut]
  s1 <- gl_sorted(clr1); s2 <- gl_sorted(clr2)
  Tc1 <- gl_top(clr1, s1, positive_only = FALSE); Tc2 <- gl_top(clr2, s2, positive_only = FALSE); rm(clr1, clr2, s1, s2)
  surrogate <- function(P) { pi <- gl_stationary(P); gl_variants(matrix(pi, G, G, byrow = TRUE), U) }
  rows <- list()
  for (clipS in c(TRUE, FALSE)) {
    sel <- Filter(function(o) o$clipS == clipS, ops); if (!length(sel)) next
    S1 <- gl_pair_score(f1, clipS); S2 <- gl_pair_score(f2, clipS)
    for (o in sel) {
      t0 <- Sys.time()
      P1 <- gl_make_P(S1, o$keff, o$cap); P2 <- gl_make_P(S2, o$keff, o$cap)
      sur1 <- surrogate(P1); sur2 <- surrogate(P2); null <- list()
      for (vn in names(sur1)) { a <- sur1[[vn]]; b <- sur2[[vn]]; sa <- gl_sorted(a); sb <- gl_sorted(b)
        null[[vn]] <- c(gl_ap_T(a, sa, gl_top(b, sb)), gl_ap_T(b, sb, gl_top(a, sa))) }   # direction 1, direction 2
      rm(sur1, sur2)
      st2 <- list(); c3_2 <- list()
      gl_walk(P2, o$depths, function(t, Pt) {                # half 2's scores, and C3 of direction 2 (P2^t against half 1's CLR targets)
        st2[[as.character(t)]] <<- gl_variants(Pt, U); c3_2[[as.character(t)]] <<- gl_c3(Pt, U$ri[Tc1], U$cj[Tc1]) })
      gl_walk(P1, o$depths, function(t, Pt) {
        v1 <- gl_variants(Pt, U); v2 <- st2[[as.character(t)]]; st2[[as.character(t)]] <<- NULL
        c3 <- list(gl_c3(Pt, U$ri[Tc2], U$cj[Tc2]), c3_2[[as.character(t)]])
        for (vn in names(v1)) {
          x1 <- v1[[vn]]; x2 <- v2[[vn]]; sv1 <- gl_sorted(x1); sv2 <- gl_sorted(x2)
          T1 <- gl_top(x1, sv1); T2 <- gl_top(x2, sv2)
          cc <- list(c(C1 = gl_ap_T(x1, sv1, Tc2), C2 = gl_ap_T(x1, sv1, T2)),    # direction 1: half 1 -> half 2
                     c(C1 = gl_ap_T(x2, sv2, Tc1), C2 = gl_ap_T(x2, sv2, T1)))   # direction 2: half 2 -> half 1
          for (d in 1:2) { ne <- n_eff[3L - d]; c3d <- c3[[d]]
            rows[[length(rows) + 1L]] <<- data.frame(direction = d, method = gl_variant_names[[vn]], keff = o$keff, cap = o$cap,
              clipS = o$clipS, clipO = gl_variant_clipO[[vn]], t = t, C1 = cc[[d]][["C1"]], C2 = cc[[d]][["C2"]], C2_null = null[[vn]][d],
              C2_adj = cc[[d]][["C2"]] - null[[vn]][d], C3 = c3d[["C3"]], C3_e10 = c3d[["C3_e10"]],
              C3_bic = ne * c3d[["C3"]] - 0.5 * GL$bic_k * log(ne), k_reach = c3d[["k_reach"]], c3_unreach = c3d[["c3_unreach"]], n_eff = ne, stringsAsFactors = FALSE) }
          rm(x1, x2, sv1, sv2)
        }
        rm(v1, v2); invisible(gc(FALSE, full = FALSE)) })
      rm(P1, P2, st2, c3_2); invisible(gc(verbose = FALSE))
      say(sprintf("operator keff %g cap %d clipS %s: %.0f s", o$keff, o$cap, o$clipS, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    rm(S1, S2)
  }
  do.call(rbind, rows)
}

## ---- operator lists ----------------------------------------------------------------------------------------------------------
# the 16 grid operators (cap = 4 x keff, at most 320) with all depths, plus the reference operator (keff 10, cap 50, clipped) at t = 9 only
gl_cap <- function(keff) as.integer(pmin(GL$cap_mult * keff, GL$cap_max))
gl_grid_ops <- function(reference = TRUE, keffs = GL$keffs, depths = GL$depths, clipS = c(TRUE, FALSE)) {
  ops <- list()
  for (cs in clipS) for (k in keffs) ops[[length(ops) + 1L]] <- list(keff = k, cap = gl_cap(k), clipS = cs, depths = depths)
  if (reference) ops[[length(ops) + 1L]] <- list(keff = GL$ref$keff, cap = GL$ref$cap, clipS = TRUE, depths = GL$ref$t)
  ops
}

## ---- I/O ---------------------------------------------------------------------------------------------------------------------
gl_read_matrix <- function(path) {                 # first column the gene name, then profiles (the expression.tsv files of the project)
  ex <- utils::read.delim(path, check.names = FALSE); X <- as.matrix(ex[, -1]); rownames(X) <- ex[[1]]; storage.mode(X) <- "double"; X
}
gl_wcsv <- function(x, path) { tmp <- paste0(path, ".tmp", Sys.getpid()); utils::write.csv(x, tmp, row.names = FALSE); file.rename(tmp, path); invisible(path) }
gl_say <- function(log = NULL) function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); if (!is.null(log)) cat(m, "\n", file = log, append = TRUE) }
gl_options <- function(args, defaults) {
  opt <- defaults; i <- 1L
  while (i <= length(args)) {
    k <- sub("^--", "", args[i]); if (!k %in% names(opt)) stop("unknown option --", k)
    v <- args[i + 1L]; opt[[k]] <- if (is.character(defaults[[k]])) v else if (is.integer(defaults[[k]])) as.integer(v) else as.numeric(v); i <- i + 2L
  }
  opt
}
gl_vmhwm_gb <- function() { l <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  x <- grep("^VmHWM:", l, value = TRUE); if (length(x)) as.numeric(gsub("[^0-9]", "", x)) / 1024^2 else NA_real_ }
