# ClrAttention — R6 orchestration over the functional core.
#
# Follows Hadley's R6 good-form rules (adv-r.hadley.nz/r6):
# UpperCamelCase class, snake_case methods, validation in $initialize(),
# informative $print(), side-effect methods return invisible(self) for
# chaining, cached stages live in private fields, read-only state is exposed
# through active bindings that fail loudly when the stage has not run.

#' Iterative CLR as attention
#'
#' Pipeline: \code{$new(data)$estimate_mi()$calibrate()$build_operator()$diffuse()}.
#' Each stage caches its result privately; \code{$mi}, \code{$clr_scores},
#' \code{$operator}, \code{$embedding} are strict read-only active bindings.
#' @export
ClrAttention <- R6::R6Class("ClrAttention",
  public = list(
    #' @description Validate and store the expression matrix.
    #' @param data numeric matrix, genes x samples. No NA/NaN/Inf.
    initialize = function(data) {
      if (!is.matrix(data) || !is.numeric(data))
        stop("data must be a numeric matrix (genes x samples)")
      if (any(!is.finite(data))) stop("data must not contain NA/NaN/Inf")
      if (nrow(data) < 2) stop("need at least 2 genes")
      if (ncol(data) < 2) stop("need at least 2 samples")
      private$data_ <- data
      private$params_ <- list(
        n_genes = nrow(data), n_samples = ncol(data),
        gene_names = rownames(data)
      )
      invisible(self)
    },

    #' @description Estimate the gene-gene MI matrix with the B-spline estimator.
    #' @param bins "fd" | "scott" | "sturges" (per-gene adaptive) or a fixed integer.
    #' @param spline_order 2 or 3.
    estimate_mi = function(bins = "fd", spline_order = 3, threads = NULL) {
      private$mi_ <- bspline_mi(private$data_, bins = bins,
                                spline_order = spline_order, threads = threads)
      # The per-gene bin counts actually used are recorded so that every
      # downstream re-estimation (permutation nulls) reuses them verbatim:
      # the null must be computed with the observed data's discretization.
      private$params_$mi <- list(bins = bins, spline_order = spline_order,
                                 threads = threads,
                                 bins_used = bins_for_genes(private$data_,
                                                            bins = bins))
      # downstream stages are stale once MI is re-estimated
      private$scores_ <- private$operator_ <- private$trajectory_ <-
        private$threshold_ <- NULL
      invisible(self)
    },

    #' @description Calibrate MI to CLR scores.
    #' @param method "normal" | "rayleigh" | "kde".
    #' @param combine "euclidean" (historical) | "stouffer" (new); normal only.
    calibrate = function(method = "normal", combine = "euclidean") {
      private$.need(private$mi_, "estimate_mi()")
      private$scores_ <- clr_calibrate(private$mi_, method = method,
                                       combine = combine)
      private$params_$calibrate <- list(method = method, combine = combine)
      private$operator_ <- private$trajectory_ <- private$threshold_ <- NULL
      invisible(self)
    },

    #' @description Select the attention threshold by permutation.
    #'
    #' Independently shuffles each gene's expression vector B times, rebuilds
    #' MI + CLR scores with the fitted settings, and pools the null scores in
    #' a streaming histogram (null matrices are never kept). P-values come
    #' from the pooled empirical null -- no parametric fit. The cutoff is
    #' chosen by Tukey's higher criticism ("hc", default) or the
    #' Benjamini-Hochberg FDR procedure ("fdr", level q). Each gene then keeps
    #' however many connections survive: per-gene degree adapts, no k to tune.
    #'
    #' Uses R's RNG: call set.seed() for reproducibility.
    #' @param B number of permutation bootstraps (default 100).
    #' @param method "hc" (default) or "fdr".
    #' @param q FDR level, used only with method = "fdr".
    #' @param threads thread count for the bootstrap MI builds (NULL = default).
    select_threshold = function(B = 100, method = c("hc", "fdr"), q = 0.05,
                                threads = NULL) {
      private$.need(private$mi_, "estimate_mi()")
      private$.need(private$scores_, "calibrate()")
      method <- match.arg(method)
      B <- as.integer(B)
      if (length(B) != 1L || is.na(B) || B < 1L)
        stop("B must be a positive integer")
      if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1)
        stop("q must be in (0, 1)")
      G <- nrow(private$data_)
      if (G < 4L) stop("select_threshold() needs at least 4 genes")
      if (!identical(private$params_$calibrate$method, "normal"))
        stop('select_threshold() requires calibrate(method = "normal"): ',
             "kde scores are non-positive log-probabilities and rayleigh ",
             "scores live in [0, 1]; the streaming null assumes z-type scores")
      M <- as.numeric(G) * (G - 1L) / 2
      mi_p <- private$params_$mi
      cal_p <- private$params_$calibrate
      X <- private$data_
      ut <- upper.tri(matrix(0, G, G))

      # Streaming pooled null: fixed-bin histogram over [0, hi), plus an
      # overflow bin. Double counts: B * M can exceed 2^31.
      nbins <- 20000L
      counts <- numeric(nbins + 1L)
      n_null <- 0
      hi <- NULL
      w <- NULL
      for (b in seq_len(B)) {
        Xp <- t(apply(X, 1L, sample))
        dimnames(Xp) <- dimnames(X)
        Sb <- clr_calibrate(
          bspline_mi(Xp, bins = mi_p$bins_used,
                     spline_order = mi_p$spline_order,
                     threads = threads),
          method = cal_p$method, combine = cal_p$combine)
        v <- Sb[ut]
        if (is.null(hi)) {
          hi <- max(1e-8, 1.25 * max(v))
          w <- hi / nbins
        }
        idx <- as.integer(v / w) + 1L
        idx[idx < 1L] <- 1L
        idx[idx > nbins] <- nbins + 1L
        counts <- counts + tabulate(idx, nbins + 1L)
        n_null <- n_null + length(v)
      }
      # suf[k] = #{null in bins >= k}. A score s in bin k gets p from
      # #{null in bins >= k}: this counts the null values sharing s's bin
      # (some of which may lie below s), so the p-value is conservative by
      # at most one bin's mass. (An earlier version used bins > k, which
      # dropped same-bin nulls >= s and was anti-conservative.) Scores in the
      # overflow bin are compared against the whole overflow count.
      suf <- rev(cumsum(rev(counts)))
      pvals <- function(s) {
        idx <- as.integer(s / w) + 1L
        idx[idx < 1L] <- 1L
        idx[idx > nbins] <- nbins + 1L
        (1 + suf[idx]) / (1 + n_null)
      }

      s_obs <- private$scores_[ut]
      p_obs <- pvals(s_obs)
      ord <- order(p_obs)
      p_sorted <- p_obs[ord]
      s_by_p <- s_obs[ord]  # scores in ascending-p-value order
      tau <- if (method == "hc") .hc_cutoff(p_sorted, s_by_p, M)
             else .fdr_cutoff(p_sorted, s_by_p, M, q)
      private$threshold_ <- tau
      private$operator_ <- private$trajectory_ <- NULL
      private$params_$threshold <- list(B = B, method = method, q = q,
                                       tau = tau, n_null = n_null)
      invisible(self)
    },

    #' @description Sparsify CLR scores and row-normalize to a stochastic operator.
    #' @param k keep the top-k scores per row (ignored if tau is given).
    #' @param tau keep scores >= tau (takes precedence over k). If NULL, uses
    #'   the threshold from select_threshold() when available, else k.
    #'   tau = Inf keeps nothing (empty operator).
    #' @param alpha diffusion weight in (0, 1].
    build_operator = function(k = 50, tau = NULL, alpha = 0.5) {
      private$.need(private$scores_, "calibrate()")
      if (is.null(tau)) tau <- private$threshold_
      S <- private$scores_
      G <- nrow(S)
      if (!is.null(tau)) {
        if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau < 0)
          stop("tau must be a single non-negative number")
        A <- S * (S >= tau)
      } else {
        k <- as.integer(k)
        if (length(k) != 1L || is.na(k) || k < 1L)
          stop("k must be a positive integer")
        k <- min(k, G - 1L)
        A <- matrix(0, G, G)
        for (i in seq_len(G)) {
          row <- S[i, ]
          cutoff <- sort(row, decreasing = TRUE)[k]
          keep <- row >= cutoff & row > 0
          A[i, keep] <- row[keep]
        }
      }
      alpha <- as.numeric(alpha)
      if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha > 1)
        stop("alpha must be in (0, 1]")
      rs <- rowSums(A)
      nz <- rs > 0
      A[nz, ] <- A[nz, ] / rs[nz]
      # Genes that attend to nobody attend to themselves: every row is then a
      # probability distribution, so P = (1-a)I + aA is stochastic and an
      # isolated gene's profile is preserved. (Previously such rows stayed
      # zero and the gene decayed as (1-a)^t toward 0; tau = Inf wiped the
      # whole embedding.)
      if (any(!nz)) A[cbind(which(!nz), which(!nz))] <- 1
      if (!is.null(rownames(S))) {
        rownames(A) <- colnames(A) <- rownames(S)
      }
      private$operator_ <- A
      private$params_$operator <- list(k = if (is.null(tau)) k else NULL,
                                       tau = tau, alpha = alpha)
      private$trajectory_ <- NULL
      invisible(self)
    },

    #' @description Iterate E <- ((1-alpha) I + alpha A_hat) E, caching the trajectory.
    #' @param steps positive integer number of diffusion steps.
    #' @param standardize if TRUE (default), E^(0) is the row-standardized
    #'   data (each gene centered and scaled to unit sd), so diffusion mixes
    #'   expression *shapes*, not baseline levels or measurement scales.
    #'   FALSE diffuses the raw matrix (the pre-2026-09-23 behavior).
    #'   Note: MI is sign-blind, so anticorrelated neighbours still partially
    #'   cancel under diffusion; see PORT_NOTES.md section 10.
    diffuse = function(steps = 10, standardize = TRUE) {
      private$.need(private$operator_, "build_operator()")
      steps <- as.integer(steps)
      if (length(steps) != 1L || is.na(steps) || steps < 1L)
        stop("steps must be a positive integer")
      alpha <- private$params_$operator$alpha
      G <- nrow(private$operator_)
      P <- (1 - alpha) * diag(G) + alpha * private$operator_
      traj <- vector("list", steps + 1L)
      E0 <- private$data_
      if (isTRUE(standardize)) {
        mu <- rowMeans(E0)
        sdv <- sqrt(rowSums((E0 - mu)^2) / (ncol(E0) - 1))
        sdv[!(sdv > 0)] <- 1  # constant gene: centered only
        E0 <- (E0 - mu) / sdv
      }
      traj[[1L]] <- E0
      for (t in seq_len(steps)) traj[[t + 1L]] <- P %*% traj[[t]]
      private$trajectory_ <- traj
      private$params_$diffuse <- list(steps = steps,
                                      standardize = isTRUE(standardize))
      invisible(self)
    },

    #' @description Drop the cached diffusion trajectory (keeps MI/scores/operator).
    reset_diffusion = function() {
      private$trajectory_ <- NULL
      private$params_$diffuse <- NULL
      invisible(self)
    },

    #' @description Summarize the pipeline state.
    print = function(...) {
      p <- private$params_
      cat("<ClrAttention>", p$n_genes, "genes x", p$n_samples, "samples\n")
      cat("  MI:        ", if (!is.null(private$mi_)) "estimated" else "--", "\n")
      cat("  CLR:       ", if (!is.null(private$scores_)) "calibrated" else "--", "\n")
      cat("  operator:  ", if (!is.null(private$operator_)) "built" else "--", "\n")
      cat("  diffusion: ",
          if (!is.null(private$trajectory_))
            paste0(length(private$trajectory_) - 1L, " steps") else "--", "\n")
      invisible(self)
    }
  ),

  active = list(
    #' @field mi the G x G MI matrix (after estimate_mi()).
    mi = function() private$.need(private$mi_, "estimate_mi()"),
    #' @field clr_scores the G x G CLR score matrix (after calibrate()).
    clr_scores = function() private$.need(private$scores_, "calibrate()"),
    #' @field threshold the permutation-selected score threshold
    #'   (after select_threshold()).
    threshold = function() private$.need(private$threshold_,
                                         "select_threshold()"),
    #' @field operator the row-stochastic attention operator (after build_operator()).
    operator = function() private$.need(private$operator_, "build_operator()"),
    #' @field embedding the diffused expression matrix (after diffuse()).
    embedding = function() {
      tr <- private$.need(private$trajectory_, "diffuse()")
      tr[[length(tr)]]
    },
    #' @field trajectory the full diffusion trajectory, E^(0..T) (after diffuse()).
    trajectory = function() private$.need(private$trajectory_, "diffuse()"),
    #' @field params read-only list of data dimensions and stage parameters.
    params = function() private$params_
  ),

  private = list(
    data_ = NULL, mi_ = NULL, scores_ = NULL, operator_ = NULL,
    trajectory_ = NULL, threshold_ = NULL, params_ = NULL,
    .need = function(value, stage) {
      if (is.null(value))
        stop("not available yet: run $", stage, " first", call. = FALSE)
      value
    }
  )
)

# Tukey's higher criticism cutoff (Donoho & Jin).
# p_sorted: ascending p-values; s_by_p: corresponding scores (descending).
# Returns the score threshold: keep pairs with p <= p_(ihat),
# i.e. scores >= s_by_p[ihat].
.hc_cutoff <- function(p_sorted, s_by_p, M) {
  n <- length(p_sorted)
  K <- max(1L, floor(n / 2))
  i <- seq_len(K)
  p <- p_sorted[i]
  ok <- p < i / M
  hc <- rep(-Inf, K)
  # p > 0 always: p-values carry a (1 + #{null > s}) / (1 + n_null) correction.
  hc[ok] <- sqrt(M) * (i[ok] / M - p[ok]) / sqrt(p[ok] * (1 - p[ok]))
  ihat <- which.max(hc)
  hcstar <- hc[ihat]
  crit <- sqrt(2 * log(log(M)))
  if (!is.finite(hcstar) || hcstar < crit)
    warning("higher criticism detected no significant signal (HC* = ",
            format(hcstar, digits = 3), " < ", format(crit, digits = 3),
            "); threshold may be noise", call. = FALSE)
  s_by_p[ihat]
}

# Benjamini-Hochberg FDR cutoff at level q. Returns Inf when nothing survives.
.fdr_cutoff <- function(p_sorted, s_by_p, M, q) {
  i <- seq_along(p_sorted)
  ok <- which(p_sorted <= i * q / M)
  if (!length(ok)) return(Inf)
  s_by_p[max(ok)]
}
