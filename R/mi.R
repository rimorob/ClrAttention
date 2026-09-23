# B-spline mutual information matrix (functional core).

#' B-spline smoothed mutual information between all gene pairs
#'
#' @param data numeric matrix, genes x samples (rows = genes).
#' @param bins "fd" (default), "scott", "sturges" for per-gene adaptive bin
#'   counts, or a single integer for a fixed count on every gene
#'   (10 reproduces the historical default).
#' @param spline_order B-spline order, 2 or 3 (historical default 3).
#' @param threads OpenMP thread count for the gene-pair loop. NULL (default)
#'   uses the core default: machine cores minus 2 (floored at 1). Must be a
#'   positive integer if given.
#' @param transform "none" (default; historical) or "rank": replace each
#'   gene by its within-gene ranks (average ties) before binning, i.e.
#'   estimate MI on the empirical copula. MI is invariant to monotone
#'   transforms, so this changes only the estimator, not the estimand; it
#'   equalizes marginal shapes across genes, which removes the
#'   heavy-tail -> more-bins -> upward-MI-bias -> false-hub pathway and makes
#'   a pooled permutation null valid for every pair.
#' @return symmetric G x G MI matrix (bits). The diagonal holds each gene's
#'   self-MI and is zeroed by [clr_calibrate()].
#' @export
bspline_mi <- function(data, bins = "fd", spline_order = 3, threads = NULL,
                       transform = c("none", "rank")) {
  transform <- match.arg(transform)
  if (!is.matrix(data) || !is.numeric(data))
    stop("data must be a numeric matrix (genes x samples)")
  if (any(!is.finite(data))) stop("data must not contain NA/NaN/Inf")
  if (nrow(data) < 2 || ncol(data) < 2)
    stop("data needs at least 2 genes and 2 samples")
  spline_order <- as.integer(spline_order)
  if (length(spline_order) != 1L || is.na(spline_order) ||
      spline_order < 2L || spline_order > 3L)
    stop("spline_order must be 2 or 3")
  if (is.null(threads)) {
    n_threads <- 0L  # core default: max(1, cores - 2)
  } else {
    n_threads <- as.integer(threads)
    if (length(n_threads) != 1L || is.na(n_threads) || n_threads < 1L)
      stop("threads must be NULL or a positive integer")
  }

  rng <- apply(data, 1L, function(x) max(x) - min(x))
  if (any(rng <= 0)) {
    bad <- which(rng <= 0)
    lab <- if (!is.null(rownames(data))) rownames(data)[bad] else bad
    stop("constant gene(s) carry no information and would become spurious ",
         "MI hubs; remove them first: ",
         paste(utils::head(lab, 10), collapse = ", "),
         if (length(bad) > 10) ", ..." else "")
  }
  data <- .transform_rows(data, transform)
  bins_vec <- bins_for_genes(data, bins = bins)
  mi <- cpp_mi_matrix(data, bins_vec, spline_order, n_threads)
  if (!is.null(rownames(data))) {
    rownames(mi) <- colnames(mi) <- rownames(data)
  }
  mi
}

# Per-gene marginal transform applied before MI estimation.
.transform_rows <- function(data, transform) {
  if (transform == "none") return(data)
  out <- t(apply(data, 1L, rank, ties.method = "average"))
  dimnames(out) <- dimnames(data)
  out
}
