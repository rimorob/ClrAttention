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
#' @return symmetric G x G MI matrix (bits). The diagonal holds each gene's
#'   self-MI and is zeroed by [clr_calibrate()].
#' @export
bspline_mi <- function(data, bins = "fd", spline_order = 3, threads = NULL) {
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

  bins_vec <- bins_for_genes(data, bins = bins)
  mi <- cpp_mi_matrix(data, bins_vec, spline_order, n_threads)
  if (!is.null(rownames(data))) {
    rownames(mi) <- colnames(mi) <- rownames(data)
  }
  mi
}
