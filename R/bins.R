# Adaptive per-gene bin counts for the B-spline MI estimator.
#
# The historical MATLAB code fixed the bin count at 10 for every gene. Here
# the count is estimated from each gene's own distribution with a classical
# optimal-width rule, then clamped to [min_bins, max_bins].

#' Optimal histogram bin count for one gene's expression vector
#'
#' @param x numeric vector (one gene's samples across conditions).
#' @param rule one of "fd" (Freedman-Diaconis, default), "scott", "sturges".
#' @param min_bins,max_bins clamp range for the adaptive count.
#' @return integer bin count in [min_bins, max_bins].
#' @details Degenerate inputs (zero IQR/sd, < 2 finite values) fall back to
#'   Sturges' rule so a positive count is always returned.
#' @export
fd_bins <- function(x, rule = c("fd", "scott", "sturges"),
                    min_bins = 5, max_bins = 50) {
  rule <- match.arg(rule)
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2) stop("fd_bins() needs at least 2 finite values")
  min_bins <- as.integer(min_bins); max_bins <- as.integer(max_bins)

  nb <- switch(rule,
    fd = {
      h <- 2 * stats::IQR(x) * n^(-1/3)
      if (!is.finite(h) || h <= 0) NA_real_ else ceiling(diff(range(x)) / h)
    },
    scott = {
      h <- 3.5 * stats::sd(x) * n^(-1/3)
      if (!is.finite(h) || h <= 0) NA_real_ else ceiling(diff(range(x)) / h)
    },
    sturges = ceiling(log2(n)) + 1
  )
  if (!is.finite(nb) || nb < 1) nb <- ceiling(log2(n)) + 1  # degenerate fallback
  as.integer(min(max(nb, min_bins), max_bins))
}

#' Per-gene bin counts for an expression matrix
#'
#' @param data numeric matrix, genes x samples.
#' @param bins "fd" | "scott" | "sturges" (per-gene adaptive), a single
#'   integer >= 2 used for every gene (historical parity mode), or an integer
#'   vector with one count per gene (used verbatim).
#' @param min_bins,max_bins clamp range applied to adaptive counts only.
#' @return integer vector of length nrow(data).
#' @export
bins_for_genes <- function(data, bins = "fd", min_bins = 5, max_bins = 50) {
  G <- nrow(data)
  if (is.character(bins)) {
    rule <- match.arg(bins, c("fd", "scott", "sturges"))
    vapply(seq_len(G), function(i) fd_bins(data[i, ], rule = rule,
                                           min_bins = min_bins,
                                           max_bins = max_bins),
           integer(1))
  } else {
    nb <- as.integer(bins)
    if (length(nb) == G && G > 1L) {
      # explicit per-gene vector (e.g. counts fitted on the observed data and
      # propagated verbatim to permutation replicates)
      if (anyNA(nb) || any(nb < 2L))
        stop("per-gene bin counts must all be integers >= 2")
      return(nb)
    }
    if (length(nb) != 1L || is.na(nb) || nb < 2L)
      stop("integer bins must be a single value >= 2 or one value per gene")
    rep(nb, G)
  }
}
