# CLR calibration of an MI matrix (functional core).
#
# method = "normal" reproduces clr.m exactly: row-wise z-scores of the
# diagonal-zeroed MI matrix (sample sd, n - 1), negatives clipped to 0 BEFORE
# bilateral combination, Euclidean combination, diagonal zeroed.
# The C++ core implements this path; rayleigh/kde are pure-R ports of clr.m.

#' Calibrate an MI matrix to CLR scores
#'
#' @param mi symmetric numeric MI matrix (e.g. from [bspline_mi()]).
#' @param method "normal" (z-scores, default), "rayleigh", or "kde".
#' @param combine bilateral combination for method = "normal": "stouffer"
#'   (default since 2026-09-23: (z_ij+z_ji)/sqrt(2), standard-normal null) or
#'   "euclidean" (historical 2007: sqrt(z_ij^2 + z_ji^2), Rayleigh-type null).
#'   Ignored (with an error if non-default) for other methods.
#' @return symmetric G x G CLR score matrix with zero diagonal.
#' @export
clr_calibrate <- function(mi, method = c("normal", "rayleigh", "kde"),
                          combine = c("stouffer", "euclidean")) {
  combine_given <- !missing(combine)
  method <- match.arg(method)
  combine <- match.arg(combine)
  if (!is.matrix(mi) || !is.numeric(mi) || nrow(mi) != ncol(mi) ||
      nrow(mi) < 2) stop("mi must be a square numeric matrix with >= 2 rows")
  if (any(!is.finite(mi))) stop("mi must not contain NA/NaN/Inf")
  if (method != "normal" && combine_given)
    stop('combine is only meaningful with method = "normal"')

  if (method == "normal") {
    code <- if (combine == "euclidean") 0L else 1L
    out <- cpp_clr_calibrate(mi, code)
  } else if (method == "rayleigh") {
    out <- .rayleigh_scores(mi)
  } else {
    out <- .kde_scores(mi)
  }
  if (!is.null(rownames(mi))) {
    rownames(out) <- colnames(out) <- rownames(mi)
  }
  out
}

# Port of clr.m 'rayleigh': per-row Rayleigh MLE (raylfit), CDF (raylcdf),
# combination A + A' - A*A'.
.rayleigh_scores <- function(mi) {
  G <- nrow(mi)
  A <- matrix(0, G, G)
  for (i in seq_len(G)) {
    sigma <- sqrt(mean(mi[i, ]^2) / 2)          # raylfit MLE
    if (!is.finite(sigma) || sigma <= 0) {
      A[i, ] <- 0
    } else {
      A[i, ] <- 1 - exp(-mi[i, ]^2 / (2 * sigma^2))  # raylcdf
    }
  }
  A <- A + t(A) - A * t(A)
  diag(A) <- 0
  A
}

# Port of clr.m 'kde': per-row Gaussian KDE on a 1000-point grid between the
# row min/max, empirical CDF from the normalized density, combination
# log(A) + log(A') with -Inf -> 0 and positives -> 0.
#
# Structural port, not bit-exact: R's density() and MATLAB's ksdensity use
# different default bandwidth selectors.
.kde_scores <- function(mi) {
  G <- nrow(mi)
  A <- matrix(0, G, G)
  for (i in seq_len(G)) {
    lo <- min(mi[i, ]); hi <- max(mi[i, ])
    if (!is.finite(lo) || hi <= lo) next
    d <- stats::density(mi[i, ], from = lo, to = hi, n = 1000)
    p <- d$y / sum(d$y)
    # sum density up to the first grid point >= v (MATLAB find(...,1,'first'))
    A[i, ] <- vapply(mi[i, ], function(v) {
      idx <- which(d$x >= v)[1L]
      if (is.na(idx)) 1 else sum(p[seq_len(idx)])
    }, numeric(1))
  }
  A <- log(A) + log(t(A))
  A[A == -Inf] <- 0
  A[A > 0] <- 0
  diag(A) <- 0
  A
}
