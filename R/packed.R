# Packed (tiled) MI -> CLR for large gene sets.
#
# clr_packed() computes B-spline MI band by band into ONE packed float32
# upper triangle held in C++ (G(G-1)/2 values: 0.65 GB at G = 18,000, against
# 2.6 GB for one dense double matrix), then calibrates it to CLR in place.
# R never holds a G x G matrix unless it asks for one: read each gene's top-k
# scores (the softmax operator needs only these), a block of genes, or the
# dense matrix. Same estimator and calibration as bspline_mi() +
# clr_calibrate(method = "normal"); the double-precision version of the same
# path is bit-identical to them (tests/testthat/test-packed.R), the float32
# scores differ by < 1e-5.

#' Packed CLR scores for large gene sets
#'
#' @param data numeric matrix, genes x samples (as for [bspline_mi()]).
#' @param bins,spline_order,threads,transform as for [bspline_mi()].
#' @param combine "stouffer" (default) or "euclidean" bilateral combination.
#' @return an object of class `clr_packed` (an external pointer to the packed
#'   scores plus the gene names); read it with [clr_packed_topk()],
#'   [clr_packed_block()], [softmax_operator()] or `as.matrix()`. The pointer
#'   does not survive saveRDS()/serialization.
#' @export
clr_packed <- function(data, bins = "hg", spline_order = 3, threads = NULL,
                       combine = c("stouffer", "euclidean"),
                       transform = c("none", "rank")) {
  combine <- match.arg(combine)
  transform <- match.arg(transform)
  if (!is.matrix(data) || !is.numeric(data))
    stop("data must be a numeric matrix (genes x samples)")
  if (any(!is.finite(data))) stop("data must not contain NA/NaN/Inf")
  if (nrow(data) < 2 || ncol(data) < 2) stop("data needs at least 2 genes and 2 samples")
  spline_order <- as.integer(spline_order)
  if (length(spline_order) != 1L || is.na(spline_order) || spline_order < 2L || spline_order > 3L)
    stop("spline_order must be 2 or 3")
  n_threads <- if (is.null(threads)) 0L else as.integer(threads)
  if (length(n_threads) != 1L || is.na(n_threads) || n_threads < 0L)
    stop("threads must be NULL or a positive integer")
  rng <- apply(data, 1L, function(x) max(x) - min(x))
  if (any(rng <= 0))
    stop("constant gene(s) carry no information; remove them first")
  data <- .transform_rows(data, transform)
  bins_vec <- bins_for_genes(data, bins = bins)
  ptr <- cpp_clr_packed(data, bins_vec, spline_order,
                        if (combine == "euclidean") 0L else 1L, n_threads)
  genes <- if (is.null(rownames(data))) as.character(seq_len(nrow(data))) else rownames(data)
  structure(list(ptr = ptr, genes = genes, n = nrow(data), combine = combine,
                 bins = bins_vec), class = "clr_packed")
}

#' @export
print.clr_packed <- function(x, ...) {
  cat(sprintf("<clr_packed> %d genes, %s combination, %.2f GB packed float32 scores\n",
              x$n, x$combine, x$n * (x$n - 1) / 2 * 4 / 1e9))
  invisible(x)
}

#' @export
as.matrix.clr_packed <- function(x, ...) {
  m <- cpp_packed_dense(x$ptr)
  dimnames(m) <- list(x$genes, x$genes)
  m
}

#' Block of a packed CLR score matrix
#' @param x a `clr_packed` object.
#' @param genes gene names or 1-based indices.
#' @return dense |genes| x |genes| matrix (diagonal 0).
#' @export
clr_packed_block <- function(x, genes) {
  idx <- if (is.character(genes)) match(genes, x$genes) else as.integer(genes)
  if (anyNA(idx)) stop("unknown genes in block")
  m <- cpp_packed_block(x$ptr, idx)
  dimnames(m) <- list(x$genes[idx], x$genes[idx])
  m
}

#' Each gene's k largest positive CLR scores
#' @param x a `clr_packed` object.
#' @param k number of partners per gene.
#' @return list(idx = G x k integer matrix of 1-based partner indices, NA
#'   where a gene has fewer than k positive scores; val = G x k scores).
#'   Ties are ordered by index, as in order(decreasing = TRUE).
#' @export
clr_packed_topk <- function(x, k = 50L) cpp_packed_topk(x$ptr, as.integer(k))

#' Softmax attention operator from packed CLR scores
#'
#' The construction of the dense path (`build_operator(softmax_keff, softmax_cap)`):
#' each gene attends over its `cap` largest positive scores with a per-row
#' temperature set by bisection so that exp(entropy) = keff; rows with at most
#' keff candidates get uniform weights; rows without candidates stay empty.
#' @param x a `clr_packed` object.
#' @param keff effective number of attended genes.
#' @param cap number of candidates per row.
#' @return sparse G x G row-stochastic matrix (dgCMatrix).
#' @export
softmax_operator <- function(x, keff = 10, cap = 50L) {
  if (!is.numeric(keff) || length(keff) != 1L || !(keff > 1)) stop("keff must be > 1")
  tk <- clr_packed_topk(x, as.integer(cap))
  G <- x$n; ii <- integer(0); jj <- integer(0); ww <- numeric(0)
  for (i in seq_len(G)) {
    ok <- !is.na(tk$idx[i, ])
    if (!any(ok)) next
    top <- tk$idx[i, ok]; s <- tk$val[i, ok]
    w <- if (length(top) <= keff) rep(1 / length(top), length(top)) else .softmax_weights(s, keff)
    ii <- c(ii, rep(i, length(top))); jj <- c(jj, top); ww <- c(ww, w)
  }
  Matrix::sparseMatrix(i = ii, j = jj, x = ww, dims = c(G, G),
                       dimnames = list(x$genes, x$genes))
}

# Per-row softmax weights with temperature by bisection (as .softmax_rows).
.softmax_weights <- function(s, keff) {
  neff <- function(lt) {
    w <- exp((s - max(s)) / exp(lt)); w <- w / sum(w)
    exp(-sum(w[w > 0] * log(w[w > 0])))
  }
  lo <- log(1e-4); hi <- log(1e4)
  for (it in 1:60) {
    mid <- (lo + hi) / 2
    if (neff(mid) < keff) lo <- mid else hi <- mid
  }
  w <- exp((s - max(s)) / exp((lo + hi) / 2))
  w / sum(w)
}
