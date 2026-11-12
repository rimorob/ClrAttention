# EXPLORATORY (docs/plans/explore_idea3_plan.md). Helpers for adaptive depth by hub decile on a log grid.
# Sourced by idea3_knee_null.R, idea3_select.R, idea3_score.R. Not for the paper or the pre-read.

I3_GRID <- c(1L, 2L, 3L, 4L, 6L, 8L, 12L, 16L, 24L, 32L, 48L, 64L, 96L, 128L)
I3_PPR  <- c(0.5, 0.7, 0.8, 0.85, 0.9, 0.95)
I3_HK   <- c(2, 4, 6, 8, 12, 16, 24, 32)
I3_TMAX <- 128L

i3_say <- function(logf, ...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

i3_source_all <- function() {
  suppressPackageStartupMessages({ library(clr); library(Matrix) })
  source("analysis/robust_z/robust_z_lib.R")       # rz_* (V0 path = package, asserted by robust_z.R)
  source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
}

# Expression matrix and 'lab' vector exactly as robust_z.R
i3_load <- function(dataset, quick = 0L) {
  if (dataset == "m3d") {
    m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
    raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
    X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], ""); rm(raw)
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
    desc <- utils::read.delim(file.path(m3d, paste0(base, ".experiment_descriptions")), stringsAsFactors = FALSE)
    exper <- sub("^.*experimenter:([^,]*).*$", "\\1", desc$description); exper[!grepl("experimenter:", desc$description)] <- NA
    lab <- tapply(exper, desc$experiment_name, function(v) v[!is.na(v)][1])[colnames(X)]
    lab[is.na(lab)] <- paste0("unknown_", which(is.na(lab)))
  } else {
    raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
    X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
    X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
    lab <- colnames(X)
  }
  storage.mode(X) <- "double"
  if (quick > 0L) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(quick)]), ]; rm(v) }
  list(X = X, lab = lab, genes = rownames(X))
}

# balanced random lab split (lab_splits.R / robust_z.R)
i3_lab_split <- function(lab) {
  sz <- table(lab); sz <- sz[sample(length(sz))]
  half <- stats::setNames(integer(length(sz)), names(sz)); tot <- c(0, 0)
  for (l in names(sz)) { h <- which.min(tot); half[l] <- h; tot[h] <- tot[h] + sz[[l]] }
  unname(half[lab])
}

i3_decile <- function(H) cut(rank(H, ties.method = "first"), 10, labels = FALSE)

# lazy-walk powers. Calls fn(Pt, t) for t = 1..tmax (Pt a dense matrix).
i3_walk <- function(P, tmax, fn) {
  Pt <- NULL
  for (t in seq_len(tmax)) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    fn(Pt, t)
    if (t %% 16L == 0L) invisible(gc(verbose = FALSE))
  }
  invisible(NULL)
}

# ---- directed degree-preserving edge swaps of the attention matrix A -------------------------
i3_rewire_cpp <- local({
  fn <- NULL
  function() {
    if (is.null(fn)) fn <<- Rcpp::cppFunction('
      IntegerVector rewire_dst(IntegerVector src, IntegerVector dst, int n, double nswap, int seed) {
        int E = src.size(); std::vector<int> d(dst.begin(), dst.end());
        std::unordered_set<long long> ex; ex.reserve((size_t)E * 2);
        for (int e = 0; e < E; e++) ex.insert((long long)src[e] * n + d[e]);
        std::mt19937_64 rng(seed); std::uniform_int_distribution<int> U(0, E - 1);
        for (double it = 0; it < nswap; it++) {
          int a = U(rng), b = U(rng);
          if (a == b) continue;
          int s1 = src[a], s2 = src[b], d1 = d[a], d2 = d[b];
          if (s1 == s2 || d1 == d2 || s1 == d2 || s2 == d1) continue;
          long long k1 = (long long)s1 * n + d2, k2 = (long long)s2 * n + d1;
          if (ex.count(k1) || ex.count(k2)) continue;
          ex.erase((long long)s1 * n + d1); ex.erase((long long)s2 * n + d2);
          ex.insert(k1); ex.insert(k2); d[a] = d2; d[b] = d1;
        }
        return wrap(d);
      }', includes = c("#include <unordered_set>", "#include <random>", "#include <vector>"))
    fn
  }
})
# A: dense G x G attention (rows sum to 1; self-loops only for empty-softmax rows, kept fixed).
i3_rewire <- function(A, seed, mult = 10) {
  G <- nrow(A); diag_e <- which(diag(A) > 0)
  B <- A; if (length(diag_e)) B[cbind(diag_e, diag_e)] <- 0
  w <- which(B != 0)                                  # column-major linear indices
  src <- (w - 1L) %% G + 1L; dst <- (w - 1L) %/% G + 1L; val <- B[w]
  nd <- i3_rewire_cpp()(as.integer(src - 1L), as.integer(dst - 1L), as.integer(G), mult * length(w), as.integer(seed)) + 1L
  R <- matrix(0, G, G); R[cbind(src, nd)] <- val
  if (length(diag_e)) R[cbind(diag_e, diag_e)] <- A[cbind(diag_e, diag_e)]
  R
}

# ---- knee of log k_i(t) (Kneedle on the grid, per gene) -----------------------------------------
# Accumulate: call i3_knee_collect(Pt, t, state) inside the walk; state$logk is G x length(grid).
i3_knee_new <- function(G) list(logk = matrix(NA_real_, G, length(I3_GRID)))
i3_knee_add <- function(st, Pt, t) { m <- match(t, I3_GRID); if (!is.na(m)) st$logk[, m] <- -log(rowSums(Pt^2)); st }
# per-gene knee: index into I3_GRID (ties: smaller t)
i3_knee_idx <- function(logk) {
  x <- log2(I3_GRID); xn <- (x - min(x)) / (max(x) - min(x))
  apply(logk, 1, function(y) { r <- range(y); if (!(r[2] > r[1])) return(length(I3_GRID))
    yn <- (y - r[1]) / (r[2] - r[1]); which.max(yn - xn) })   # which.max returns the first maximum
}
# snap a log2 depth to the nearest grid depth (ties: smaller t)
i3_snap <- function(l2) { d <- abs(log2(I3_GRID) - l2); I3_GRID[which.min(round(d, 12))] }

# ---- selection rule on a mean criterion curve ----------------------------------------------------
# crit: numeric over `pos` (log2 positions of the candidates, increasing); returns list(idx, lo, hi, width, argmax)
i3_pick <- function(crit, pos) {
  k <- which.max(crit); thr <- 0.99 * crit[k]
  lo <- k; while (lo > 1L && crit[lo - 1L] >= thr) lo <- lo - 1L
  hi <- k; while (hi < length(crit) && crit[hi + 1L] >= thr) hi <- hi + 1L
  if (hi == lo) idx <- k else { ctr <- (pos[lo] + pos[hi]) / 2; d <- round(abs(pos - ctr), 12); idx <- which.min(d) }   # equal distance: smaller
  list(idx = idx, argmax = k, lo = lo, hi = hi, width_oct = pos[hi] - pos[lo])
}

# ---- walk-mass readouts --------------------------------------------------------------------------
# R: row i of P^{t_i}; S = (R + R')/2, zero diagonal.  Pt_fun(t) returns the dense P^t.
i3_symmetrise <- function(R) { S <- (R + t(R)) / 2; diag(S) <- 0; S }
