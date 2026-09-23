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
    #' @param threads OpenMP threads (NULL = cores - 2).
    #' @param transform "none" (default; MI on raw values) or "rank" (empirical copula);
    #'   see [bspline_mi()].
    estimate_mi = function(bins = "hg", spline_order = 3,
                           threads = NULL, transform = "none") {
      private$mi_ <- bspline_mi(private$data_, bins = bins,
                                spline_order = spline_order, threads = threads,
                                transform = transform)
      # The per-gene bin counts actually used are recorded so that every
      # downstream re-estimation (permutation nulls) reuses them verbatim:
      # the null must be computed with the observed data's discretization.
      private$params_$mi <- list(bins = bins, spline_order = spline_order,
                                 threads = threads, transform = transform,
                                 bins_used = bins_for_genes(
                                   .transform_rows(private$data_, transform),
                                   bins = bins))
      # downstream stages are stale once MI is re-estimated
      private$scores_ <- private$operator_ <- private$trajectory_ <-
        private$threshold_ <- NULL
      invisible(self)
    },

    #' @description Calibrate MI to CLR scores.
    #' @param method "normal" | "rayleigh" | "kde".
    #' @param combine "stouffer" (default) | "euclidean" (historical 2007);
    #'   normal only.
    calibrate = function(method = "normal", combine = "stouffer") {
      private$.need(private$mi_, "estimate_mi()")
      private$scores_ <- if (method == "normal")
        clr_calibrate(private$mi_, method = method, combine = combine) else
        clr_calibrate(private$mi_, method = method)
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
    #' chosen by the Benjamini-Hochberg FDR procedure ("fdr", level q,
    #' default) or permutation-calibrated higher criticism ("hc"). Each gene
    #' then keeps
    #' however many connections survive: per-gene degree adapts, no k to tune.
    #'
    #' Uses R's RNG: call set.seed() for reproducibility.
    #' @param B number of permutation bootstraps (default 100).
    #' @param method "fdr" (default; Benjamini-Hochberg at level q) or "hc"
    #'   (permutation-calibrated higher criticism).
    #' @param q FDR level, used only with method = "fdr".
    #' @param threads thread count for the bootstrap MI builds (NULL = default).
    #' @param hc_alpha0 HC search range: only the top hc_alpha0 * M pairs are
    #'   candidate cutoffs (Donoho & Jin 2004/2008 use alpha0 <= 1/2; small
    #'   values stop the argmax from wandering into the null bulk).
    #' @param hc_level significance level of the permutation-calibrated HC
    #'   gate: the observed HC* must exceed the (1 - hc_level) quantile of the
    #'   leave-one-out HC* of the B permutation replicates, otherwise no edge
    #'   is selected (tau = Inf). The asymptotic sqrt(2 log log M) bound is
    #'   NOT used as a gate: it is the typical size of null HC*, not a
    #'   significance cutoff (~25% of null HC* exceed it at M ~ 2000).
    #' @param statistic which statistic the permutation null is built for.
    #'   "clr" (default): edges are selected on their CLR score, i.e. on
    #'   being exceptional *relative to each gene's own background*. This is
    #'   what makes the selection robust to global dependence (growth-rate or
    #'   batch programs that touch every gene): on a synthetic with a global
    #'   factor (loading 0.35), BH on the CLR null kept 1.2% of pairs at
    #'   precision 0.99, while BH on the MI null kept 45% at precision 0.03.
    #'   Caveat: CLR scores are not exactly pivotal between observed and
    #'   permuted data. When modules are large relative to G (toy data), a
    #'   gene's own module inflates its row background and the CLR null
    #'   becomes conservative, possibly selecting nothing; this is a small-G
    #'   artifact, not the regime CLR was designed for.
    #'   "mi": the null is built on raw MI (CLR remains the attention weight;
    #'   a selected pair whose clipped CLR score is 0 gets zero weight, so
    #'   its genes can still end up isolated, i.e. self-loop only).
    #'   Exact per pair under the rank transform with equal bins, but it tests
    #'   "any dependence", which in real compendia is true of most pairs;
    #'   use only for small or confounder-free data.
    select_threshold = function(B = 100, method = c("fdr", "hc"), q = 0.05,
                                threads = NULL, hc_alpha0 = 0.1,
                                hc_level = 0.05,
                                statistic = c("clr", "mi")) {
      statistic <- match.arg(statistic)
      private$.need(private$mi_, "estimate_mi()")
      private$.need(private$scores_, "calibrate()")
      method <- match.arg(method)
      B <- as.integer(B)
      if (length(B) != 1L || is.na(B) || B < 1L)
        stop("B must be a positive integer")
      if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1)
        stop("q must be in (0, 1)")
      if (!is.numeric(hc_alpha0) || length(hc_alpha0) != 1L ||
          !(hc_alpha0 > 0 && hc_alpha0 <= 0.5))
        stop("hc_alpha0 must be in (0, 0.5]")
      if (!is.numeric(hc_level) || length(hc_level) != 1L ||
          !(hc_level > 0 && hc_level < 1))
        stop("hc_level must be in (0, 1)")
      if (method == "hc" && B < 2L)
        stop("method = \"hc\" needs B >= 2 (leave-one-out HC calibration)")
      if (method == "hc" && (B + 1) * hc_level < 1)
        warning("B = ", B, " is too small to resolve hc_level = ", hc_level,
                "; use B >= ", ceiling(1 / hc_level - 1), call. = FALSE)
      G <- nrow(private$data_)
      if (G < 4L) stop("select_threshold() needs at least 4 genes")
      if (statistic == "mi" && !identical(private$params_$mi$transform, "rank"))
        warning("statistic = \"mi\" without transform = \"rank\": the null ",
                "distribution of B-spline MI depends on the two marginal ",
                "shapes, so a pooled MI null is not exact per pair (skewed ",
                "or heavy-tailed genes get anti-conservative p-values)",
                call. = FALSE)
      if (statistic == "mi" && length(unique(private$params_$mi$bins_used)) > 1L)
        warning("statistic = \"mi\" with unequal per-gene bin counts: the ",
                "pooled MI null is not exact per pair (MI bias depends on ",
                "bins_i x bins_j); use transform = \"rank\" or a fixed bin ",
                "count", call. = FALSE)
      if (statistic == "clr" &&
          !identical(private$params_$calibrate$method, "normal"))
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
      rep_counts <- matrix(0, nbins + 1L, B)  # per-replicate histograms (HC gate)
      n_null <- 0
      hi <- NULL
      w <- NULL
      for (b in seq_len(B)) {
        Xp <- t(apply(X, 1L, sample))
        dimnames(Xp) <- dimnames(X)
        Mb <- bspline_mi(Xp, bins = mi_p$bins_used,
                         spline_order = mi_p$spline_order,
                         threads = threads, transform = mi_p$transform)
        v <- if (statistic == "mi") Mb[ut] else
          clr_calibrate(Mb, method = cal_p$method,
                        combine = cal_p$combine)[ut]
        if (is.null(hi)) {
          hi <- max(1e-8, 1.25 * max(v))
          w <- hi / nbins
        }
        idx <- as.integer(pmin(v / w, nbins)) + 1L
        idx[idx < 1L] <- 1L
        idx[idx > nbins] <- nbins + 1L
        cb <- tabulate(idx, nbins + 1L)
        rep_counts[, b] <- cb
        counts <- counts + cb
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
        idx <- as.integer(pmin(s / w, nbins)) + 1L
        idx[idx < 1L] <- 1L
        idx[idx > nbins] <- nbins + 1L
        (1 + suf[idx]) / (1 + n_null)
      }

      # TODO(empirical null, Efron 2004 JASA 99:96-104): implement an
      # alternative to the full-independence permutation null for the CLR
      # statistic -- estimate the null (center, scale) from the central bulk
      # of the OBSERVED CLR z-scores (e.g. Efron's central-matching /
      # locfdr-style fit) and compute p-values or local fdr against it. This
      # follows CLR's own logic (the data are their own background) and is
      # the principled fix if the permutation CLR null proves miscalibrated
      # at compendium scale (see CITATION_LOG.md D19; small-G compression).
      # Requested by the user 2026-09-23; not yet implemented.
      s_obs <- if (statistic == "mi") private$mi_[ut] else private$scores_[ut]
      hc_info <- NULL
      if (method == "hc") {
        idx_obs <- as.integer(pmin(s_obs / w, nbins)) + 1L
        idx_obs[idx_obs < 1L] <- 1L
        idx_obs[idx_obs > nbins] <- nbins + 1L
        obs_counts <- tabulate(idx_obs, nbins + 1L)
        obs <- .hc_binned(obs_counts, suf, n_null, M, hc_alpha0)
        # Leave-one-out null HC*: replicate b scored against the pool of the
        # other B - 1 replicates, exactly as the observed data are scored
        # against a pool it is not part of.
        hc_null <- vapply(seq_len(B), function(b) {
          cb <- rep_counts[, b]
          .hc_binned(cb, suf - rev(cumsum(rev(cb))), n_null - M, M,
                     hc_alpha0)$hc_star
        }, numeric(1))
        kq <- min(B, ceiling((B + 1) * (1 - hc_level)))
        crit <- sort(hc_null)[kq]
        if (is.finite(obs$hc_star) && obs$hc_star > crit) {
          # NB: the HC tau is a histogram bin's lower edge (keeps every score
          # in bins >= k), not an observed score as with BH.
          tau <- (obs$k - 1L) * w
        } else {
          tau <- Inf
          message("higher criticism: no signal beyond the permutation null ",
                  "(HC* = ", format(obs$hc_star, digits = 3),
                  " <= null ", 100 * (1 - hc_level), "% quantile ",
                  format(crit, digits = 3), "); no edges selected")
        }
        hc_info <- list(hc_star = obs$hc_star, hc_crit = crit,
                        hc_null = hc_null, alpha0 = hc_alpha0,
                        level = hc_level)
      } else {
        p_obs <- pvals(s_obs)
        ord <- order(p_obs)
        tau <- .fdr_cutoff(p_obs[ord], s_obs[ord], M, q)
      }
      private$threshold_ <- tau
      private$null_ <- list(suf = suf, w = w, n_null = n_null, nbins = nbins,
                            M = M, statistic = statistic)
      private$operator_ <- private$trajectory_ <- NULL
      private$params_$threshold <- list(B = B, method = method, q = q,
                                       statistic = statistic,
                                       tau = tau, n_null = n_null,
                                       hc = hc_info)
      invisible(self)
    },

    #' @description Sparsify CLR scores and row-normalize to a stochastic operator.
    #' @param k keep the top-k scores per row (ignored if tau is given).
    #' @param tau keep scores >= tau (takes precedence over k). If NULL, uses
    #'   the threshold from select_threshold() when available, else k.
    #'   tau = Inf keeps nothing (empty operator).
    #' @param alpha diffusion weight in (0, 1].
    build_operator = function(k = 50, tau = NULL, alpha = 0.5,
                              topk_union = NULL, softmax_keff = NULL,
                              softmax_cap = 50L) {
      private$.need(private$scores_, "calibrate()")
      S <- private$scores_
      G <- nrow(S)
      selected <- NULL
      if (!is.null(softmax_keff)) {
        A <- .softmax_rows(S, keff = softmax_keff, cap = softmax_cap)
        sel_label <- sprintf("softmax_keff%g_cap%d", softmax_keff,
                             as.integer(softmax_cap))
      } else {
        if (is.null(tau) && !is.null(private$threshold_)) {
          selected <- self$edges    # permutation-selected edge set
          A <- S * selected
          sel_label <- private$params_$threshold$statistic
        } else if (!is.null(tau)) {
          if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau < 0)
            stop("tau must be a single non-negative number")
          A <- S * (S >= tau)
          diag(A) <- 0
          sel_label <- "clr_tau"
        } else {
          A <- .topk_rows(S, k)
          sel_label <- "top_k"
        }
        if (!is.null(topk_union)) {
          # Directed: gene i additionally attends to its own top-k CLR
          # neighbours (row i of S), whether or not they passed selection.
          Tk <- .topk_rows(S, topk_union)
          A <- ifelse(Tk > 0, Tk, A)
          sel_label <- paste0(sel_label, "+top", as.integer(topk_union))
        }
      }
      private$.finish_operator(A, alpha, list(
        k = if (sel_label == "top_k") as.integer(k) else NULL,
        tau = if (!is.null(selected)) private$threshold_ else tau,
        topk_union = topk_union, softmax_keff = softmax_keff,
        selection = sel_label))
    },

    #' @description Install an externally built attention matrix (e.g. a
    #'   control operator from another similarity). Rows are normalized to
    #'   sum to 1; rows with no mass get a self-loop.
    #' @param A non-negative G x G matrix (row i = gene i's attention).
    #' @param alpha diffusion weight in (0, 1].
    #' @param label free-text description stored in params$operator.
    set_operator = function(A, alpha = 0.5, label = "external") {
      G <- nrow(private$data_)
      if (!is.matrix(A) || !is.numeric(A) || any(dim(A) != G))
        stop("A must be a numeric ", G, " x ", G, " matrix")
      if (any(!is.finite(A)) || any(A < 0)) stop("A must be finite and >= 0")
      private$.finish_operator(A, alpha, list(selection = label))
    },

    #' @description Restore a saved selection (null distribution + threshold)
    #'   without re-running permutations, e.g. from a saved run.
    #' @param null a list as returned by the null_distribution field.
    #' @param tau the saved threshold.
    #' @param q the saved FDR level.
    restore_threshold = function(null, tau, q = 0.05) {
      private$.need(private$scores_, "calibrate()")
      stopifnot(is.list(null), all(c("suf", "w", "n_null", "nbins", "M",
                                     "statistic") %in% names(null)))
      private$null_ <- null
      private$threshold_ <- tau
      private$params_$threshold <- list(method = "fdr", q = q, tau = tau,
                                        statistic = null$statistic,
                                        n_null = null$n_null, restored = TRUE)
      private$operator_ <- private$trajectory_ <- NULL
      invisible(self)
    },

    #' @description Re-run the Benjamini-Hochberg selection at a new FDR level
    #'   against the permutation null stored by select_threshold() -- no new
    #'   permutations.
    #' @param q FDR level in (0, 1).
    reselect = function(q) {
      nl <- private$.need(private$null_, "select_threshold()")
      if (!is.numeric(q) || length(q) != 1L || !(q > 0 && q < 1))
        stop("q must be in (0, 1)")
      if (!(nl$n_null > 0)) {
        if (isTRUE(all.equal(q, private$params_$threshold$q))) {
          private$operator_ <- private$trajectory_ <- NULL
          return(invisible(self))    # restored selection without a null
        }
        stop("no stored permutation null: cannot reselect at a new q")
      }
      G <- nrow(private$data_)
      ut <- upper.tri(matrix(0, G, G))
      s_obs <- if (nl$statistic == "mi") private$mi_[ut] else private$scores_[ut]
      idx <- as.integer(pmin(s_obs / nl$w, nl$nbins)) + 1L
      idx[idx < 1L] <- 1L
      p_obs <- (1 + nl$suf[idx]) / (1 + nl$n_null)
      ord <- order(p_obs)
      tau <- .fdr_cutoff(p_obs[ord], s_obs[ord], nl$M, q)
      private$threshold_ <- tau
      private$params_$threshold$method <- "fdr"
      private$params_$threshold$q <- q
      private$params_$threshold$tau <- tau
      private$operator_ <- private$trajectory_ <- NULL
      invisible(self)
    },

    #' @description Diffuse expression through the attention operator,
    #'   E_i <- (1-alpha) E_i + alpha * sum_j A_ij v_{j->i}(E_j), caching the
    #'   trajectory.
    #' @param steps positive integer number of diffusion steps.
    #' @param standardize if TRUE (default), E^(0) is the row-standardized
    #'   data, so diffusion mixes expression shapes, not levels or scales.
    #'   Required (forced) for values = "signed" / "conditional".
    #' @param values the value passed from attended gene j to gene i -- the
    #'   continuous analogue of the transformer's W_V:
    #'   "raw": v = E_j (linear diffusion; MI is sign-blind, so a repressed
    #'   target and its repressor partially cancel);
    #'   "signed": v = sign(cor(x_i, x_j)) E_j (fixes sign, still linear);
    #'   "conditional": v = f_ij(E_j), where f_ij(u) = E[z_i | z_j = u] is the
    #'   nonparametric regression of gene i on gene j read off the same
    #'   B-spline basis the MI estimator uses (gene j's bin count and spline
    #'   order): f_ij(u) = sum_b B_jb(u) mu_ijb with
    #'   mu_ijb = sum_s B_jb(z_js) z_is / sum_s B_jb(z_js). It carries sign,
    #'   non-monotone shape (e.g. quadratic dependence) and strength
    #'   (f_ij is ~0 for weakly dependent pairs), with no learned parameters.
    #'   Values are pair-specific (they depend on the query i as well as the
    #'   key j), like edge-conditioned messages in graph networks. f_ij is
    #'   fitted once on the data and re-applied to the diffused profiles at
    #'   every step (values outside z_j's range are clamped to it).
    diffuse = function(steps = 10, standardize = TRUE,
                       values = c("raw", "signed", "conditional")) {
      private$.need(private$operator_, "build_operator()")
      values <- match.arg(values)
      steps <- as.integer(steps)
      if (length(steps) != 1L || is.na(steps) || steps < 1L)
        stop("steps must be a positive integer")
      if (values != "raw") standardize <- TRUE
      alpha <- private$params_$operator$alpha
      A <- private$operator_
      G <- nrow(A)
      E0 <- private$data_
      if (isTRUE(standardize)) {
        mu <- rowMeans(E0)
        sdv <- sqrt(rowSums((E0 - mu)^2) / (ncol(E0) - 1))
        sdv[!(sdv > 0)] <- 1  # constant gene: centered only
        E0 <- (E0 - mu) / sdv
      }
      traj <- vector("list", steps + 1L)
      traj[[1L]] <- E0
      if (values == "raw") {
        P <- Matrix::Matrix((1 - alpha) * diag(G) + alpha * A, sparse = TRUE)
        for (t in seq_len(steps))
          traj[[t + 1L]] <- as.matrix(P %*% traj[[t]])
      } else if (values == "signed") {
        off <- A > 0 & row(A) != col(A)
        R <- stats::cor(t(E0))
        As <- A
        As[off] <- A[off] * sign(R[off])
        P <- Matrix::Matrix((1 - alpha) * diag(G) + alpha * As, sparse = TRUE)
        for (t in seq_len(steps))
          traj[[t + 1L]] <- as.matrix(P %*% traj[[t]])
      } else {
        vm <- private$.value_model(E0)
        for (t in seq_len(steps))
          traj[[t + 1L]] <- (1 - alpha) * traj[[t]] +
            alpha * private$.messages(traj[[t]], vm)
      }
      private$trajectory_ <- traj
      private$params_$diffuse <- list(steps = steps,
                                      standardize = isTRUE(standardize),
                                      values = values)
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
    #' @field edges logical G x G matrix of permutation-selected pairs
    #'   (symmetric, FALSE diagonal): MI >= threshold for statistic = "mi",
    #'   CLR >= threshold for statistic = "clr".
    edges = function() {
      tau <- private$.need(private$threshold_, "select_threshold()")
      st <- private$params_$threshold$statistic
      V <- if (identical(st, "mi")) private$mi_ else private$scores_
      E <- V >= tau
      diag(E) <- FALSE
      E
    },
    #' @field operator the row-stochastic attention operator (after build_operator()).
    operator = function() private$.need(private$operator_, "build_operator()"),
    #' @field embedding the diffused expression matrix (after diffuse()).
    embedding = function() {
      tr <- private$.need(private$trajectory_, "diffuse()")
      tr[[length(tr)]]
    },
    #' @field trajectory the full diffusion trajectory, E^(0..T) (after diffuse()).
    trajectory = function() private$.need(private$trajectory_, "diffuse()"),
    #' @field null_distribution the pooled permutation null stored by
    #'   select_threshold() (histogram survival counts, bin width, size,
    #'   statistic), for reuse by reselect() or restore_threshold().
    null_distribution = function() private$.need(private$null_,
                                                  "select_threshold()"),
    #' @field params read-only list of data dimensions and stage parameters.
    params = function() private$params_
  ),

  private = list(
    data_ = NULL, mi_ = NULL, scores_ = NULL, operator_ = NULL,
    trajectory_ = NULL, threshold_ = NULL, params_ = NULL, null_ = NULL,

    # Row-normalize, give empty rows a self-loop, store, invalidate trajectory.
    .finish_operator = function(A, alpha, info) {
      alpha <- as.numeric(alpha)
      if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha > 1)
        stop("alpha must be in (0, 1]")
      rs <- rowSums(A)
      nz <- rs > 0
      A[nz, ] <- A[nz, ] / rs[nz]
      # Genes that attend to nobody attend to themselves: every row is then a
      # probability distribution, so P = (1-a)I + aA is stochastic and an
      # isolated gene's profile is preserved.
      if (any(!nz)) A[cbind(which(!nz), which(!nz))] <- 1
      if (!is.null(rownames(private$data_)))
        rownames(A) <- colnames(A) <- rownames(private$data_)
      private$operator_ <- A
      private$params_$operator <- c(info, list(alpha = alpha))
      private$trajectory_ <- NULL
      invisible(self)
    },
    .need = function(value, stage) {
      if (is.null(value))
        stop("not available yet: run $", stage, " first", call. = FALSE)
      value
    },

    # Conditional-expectation value model on standardized data Z:
    # for every key gene j with out-edges (A[i, j] > 0, i != j), the basis
    # range of z_j and, per attending query i, mu_ij = E[z_i | bin of z_j].
    .value_model = function(Z) {
      A <- private$operator_
      G <- nrow(A)
      nb <- private$params_$mi$bins_used
      k <- as.integer(private$params_$mi$spline_order)
      keys <- lapply(seq_len(G), function(j) {
        q <- which(A[, j] > 0 & seq_len(G) != j)
        if (!length(q)) return(NULL)
        z <- Z[j, ]
        lo <- min(z); hi <- max(z)
        if (!(hi > lo)) return(NULL)
        W <- cpp_weights_at(z, lo, hi, k, nb[j])        # N x nb_j
        den <- colSums(W)
        num <- crossprod(W, t(Z[q, , drop = FALSE]))    # nb_j x |q|
        mu <- num / ifelse(den > 1e-12, den, Inf)       # empty bin -> 0
        list(j = j, q = q, w = A[q, j], lo = lo, hi = hi, nb = nb[j],
             mu = mu)
      })
      list(keys = keys[!vapply(keys, is.null, NA)], k = k,
           self = diag(A))
    },

    # Messages sum_j A_ij f_ij(E_j) for all i (G x N); self-loops pass E_i.
    .messages = function(E, vm) {
      M <- E * vm$self
      for (key in vm$keys) {
        W <- cpp_weights_at(E[key$j, ], key$lo, key$hi, vm$k, key$nb)
        pred <- W %*% key$mu                              # N x |q|
        M[key$q, ] <- M[key$q, ] + t(pred) * key$w
      }
      M
    }
  )
)

# Tukey's higher criticism on a binned score histogram (Donoho & Jin 2004;
# HC thresholding: Donoho & Jin 2008).
# obs_counts: histogram of the M scores under test (same bins as the null);
# suf: null survival counts, suf[k] = #{null in bins >= k}; n_null: null size.
# Candidate cutoffs are bin lower edges. Keeping bins >= k keeps
# i_k = #{scores in bins >= k} pairs, whose largest p-value is
# p_k = (1 + suf[k]) / (1 + n_null). HC(k) = sqrt(M)(i_k/M - p_k) /
# sqrt(p_k(1 - p_k)), maximized over i_k <= max(alpha0 * M, 10). The HC+
# floor p_k >= 1/M is NOT applied: with a permutation null the p-value
# resolution is 1/(B*M) << 1/M, so the floor would discard exactly the
# strongest edges (all of them, when M is small); the extreme-p blow-up it
# guards against is instead absorbed by the leave-one-out permutation
# calibration of HC*, whose null replicates blow up the same way.
# Evaluating at bin edges treats tied scores as one block.
# Returns list(hc_star, k): the maximum and the bin index attaining it.
.hc_binned <- function(obs_counts, suf, n_null, M, alpha0) {
  i_k <- rev(cumsum(rev(obs_counts)))
  p_k <- (1 + suf) / (1 + n_null)
  ok <- i_k >= 1 & i_k <= max(alpha0 * M, 10) & p_k < 1
  if (!any(ok)) return(list(hc_star = -Inf, k = NA_integer_))
  hc <- rep(-Inf, length(i_k))
  hc[ok] <- sqrt(M) * (i_k[ok] / M - p_k[ok]) / sqrt(p_k[ok] * (1 - p_k[ok]))
  k <- which.max(hc)
  list(hc_star = hc[k], k = k)
}

# Benjamini-Hochberg FDR cutoff at level q. Returns Inf when nothing survives.
.fdr_cutoff <- function(p_sorted, s_by_p, M, q) {
  i <- seq_along(p_sorted)
  ok <- which(p_sorted <= i * q / M)
  if (!length(ok)) return(Inf)
  # Histogram p-values are tied within a bin; order() breaks those ties by
  # index, so the score at position max(ok) is an arbitrary member of the
  # last rejected block. Return the smallest score among ALL rejected p's,
  # so score >= tau keeps exactly the BH rejection set.
  p_cut <- p_sorted[max(ok)]
  min(s_by_p[p_sorted <= p_cut])
}

# Directed top-k: row i keeps its k largest positive off-diagonal scores.
.topk_rows <- function(S, k) {
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L) stop("k must be a positive integer")
  G <- nrow(S)
  k <- min(k, G - 1L)
  A <- matrix(0, G, G)
  for (i in seq_len(G)) {
    row <- S[i, ]; row[i] <- -Inf
    top <- order(row, decreasing = TRUE)[seq_len(k)]
    top <- top[row[top] > 0]
    A[i, top] <- row[top]
  }
  A
}

# Softmax attention over each row's top-`cap` positive scores, with a
# per-row temperature chosen so the effective number of attended genes,
# exp(entropy), equals keff (bisection on log temperature). Rows with fewer
# than keff candidates get uniform weights over what they have.
.softmax_rows <- function(S, keff = 10, cap = 50L) {
  if (!is.numeric(keff) || length(keff) != 1L || !(keff > 1))
    stop("softmax_keff must be > 1")
  cap <- as.integer(cap)
  G <- nrow(S)
  A <- matrix(0, G, G)
  for (i in seq_len(G)) {
    row <- S[i, ]; row[i] <- -Inf
    top <- order(row, decreasing = TRUE)[seq_len(min(cap, G - 1L))]
    top <- top[row[top] > 0]
    if (!length(top)) next
    s <- row[top]
    if (length(top) <= keff) { A[i, top] <- 1 / length(top); next }
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
    A[i, top] <- w / sum(w)
  }
  A
}
