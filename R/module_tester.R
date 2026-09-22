# ModuleTester — R6 alignment of found modules against reference regulons.
#
# Downstream, the reference is RegulonDB: each regulon is the TF plus its
# target gene set, and a found module is scored by its best-matching regulon.
# For now the framework is exercised on the synthetic planted modules, where
# the reference is ground truth. The tester is deliberately ID-agnostic:
# modules are character vectors of gene identifiers in one shared ID space.
#
# Scoring starts from the pairwise set-F1 matrix between reference and found
# modules:
#   precision = |R cap F| / |F|, recall = |R cap F| / |R|,
#   F1 = 2PR / (P + R), Jaccard = |R cap F| / |R cup F|.
#
# Alignment (reference -> found) comes in three flavors, chosen by method:
#   "soft" (default): entropic-regularized assignment. A doubly stochastic
#     matrix P maximizes <P, F1> + epsilon * H(P) (Sinkhorn); each reference
#     module's reported score is the *expected* F1 under its assignment
#     distribution, sum_j P_ij F1_ij / sum_j P_ij. This is the probabilistic
#     generalization of one-to-one matching: epsilon -> 0 recovers the hard
#     Hungarian assignment, epsilon -> Inf spreads mass uniformly. It keeps
#     the one-to-one competition (a found module's mass is split across the
#     regulons that claim it) while staying smooth and differentiable.
#   "hungarian": hard one-to-one assignment maximizing total F1
#     (Kuhn-Munkres). Each found module explains at most one regulon.
#   "best": each regulon independently takes its argmax-F1 found module
#     (many-to-one allowed). Kept for comparison; gameable by hairballs.
#
# Rectangular problems are padded with dummy modules (F1 = 0) so the
# assignment is square; mass placed on dummies scores 0, i.e. a regulon no
# found module claims is a clean miss.

# Numerically stable log-sum-exp of a vector.
.logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Entropic-regularized assignment via Sinkhorn iterations, computed in the
# log domain for stability at small epsilon.
# F: n x n reward matrix (here: set F1, in [0,1]).
# tol is on the max absolute deviation of row/column sums from 1/n. The
# default 1e-5 sits comfortably below the ~1e-6 floating-point stagnation
# floor of the log-domain iteration on ill-conditioned (near-permutation)
# kernels, yet far below any statistically meaningful probability error:
# expected F1 scores are unaffected at that level, and .compute_alignment()
# row-renormalizes regardless. Genuine nonconvergence (e.g. max_iter hit
# far from a fixed point) still warns loudly.
# Returns a list: P (doubly stochastic matrix with uniform 1/n margins,
# maximizing <P, F> + epsilon * H(P)), converged, iters, marginal_error
# (max deviation of row/column sums from 1/n at termination).
.sinkhorn <- function(F, epsilon = 0.1, max_iter = 1000L, tol = 1e-5) {
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0)
    stop("epsilon must be a single finite positive number", call. = FALSE)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0)
    stop("tol must be a single finite positive number", call. = FALSE)
  if (length(max_iter) != 1L || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != floor(max_iter))
    stop("max_iter must be a single positive integer", call. = FALSE)
  max_iter <- as.integer(max_iter)
  n <- nrow(F)
  if (n < 1L) stop("F must have at least one row", call. = FALSE)
  # Shifted log-kernel: entries in (-Inf, 0]; log-sum-exp absorbs the -Inf
  # terms, so rows/columns that underflow entirely still behave sanely.
  logK <- F / epsilon
  logK <- logK - max(logK)
  logr <- -log(n)   # target log-marginal: uniform 1/n
  f <- rep(0, n)    # log row potentials
  g <- rep(0, n)    # log column potentials
  converged <- FALSE
  it <- 0L
  err <- Inf
  # logK[i, j] + g[j]: rep(g, each = n) aligns g with columns (col-major).
  # logK[i, j] + f[i]: plain f recycles down columns (col-major).
  for (it in seq_len(max_iter)) {
    f <- logr - vapply(seq_len(n),
                       function(i) .logsumexp(logK[i, ] + g), numeric(1))
    g <- logr - vapply(seq_len(n),
                       function(j) .logsumexp(logK[, j] + f), numeric(1))
    logP <- logK + f + rep(g, each = n)
    P <- exp(logP)
    err <- max(abs(rowSums(P) - 1 / n), abs(colSums(P) - 1 / n))
    if (err < tol) { converged <- TRUE; break }
  }
  if (!converged)
    warning("Sinkhorn did not converge in ", max_iter,
            " iterations (marginal error ", signif(err, 3), ")", call. = FALSE)
  list(P = P, converged = converged, iters = it, marginal_error = err)
}

# Kuhn-Munkres (Hungarian) for the linear assignment problem, O(n^3).
# cost: n x n matrix; returns integer vector a with a[j] = row assigned to
# column j, minimizing total cost. (cp-algorithms formulation, 1-indexed.)
.hungarian <- function(cost) {
  n <- nrow(cost)
  u <- numeric(n + 1); v <- numeric(n + 1)
  p <- integer(n + 1); way <- integer(n + 1)
  for (i in seq_len(n)) {
    p[1] <- i
    j0 <- 1L
    minv <- rep(Inf, n + 1)
    used <- rep(FALSE, n + 1)
    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- 1L
      for (j in 2:(n + 1)) {
        if (!used[j]) {
          cur <- cost[i0, j - 1] - u[i0] - v[j]
          if (cur < minv[j]) { minv[j] <- cur; way[j] <- j0 }
          if (minv[j] < delta) { delta <- minv[j]; j1 <- j }
        }
      }
      for (j in seq_len(n + 1)) {
        if (used[j]) { u[p[j]] <- u[p[j]] + delta; v[j] <- v[j] - delta }
        else minv[j] <- minv[j] - delta
      }
      j0 <- j1
      if (p[j0] == 0L) break
    }
    repeat {
      j1 <- way[j0]; p[j0] <- p[j1]; j0 <- j1
      if (j0 == 1L) break
    }
  }
  a <- integer(n)
  for (j in 2:(n + 1)) a[j - 1] <- p[j]
  a  # a[col] = row
}

# Pairwise set statistics between every reference and found module.
# Returns a list of n_ref x n_found matrices: f1, precision, recall, jaccard.
.pairwise_set_stats <- function(ref, found) {
  nr <- length(ref); nf <- length(found)
  F1 <- P <- R <- J <- matrix(0, nr, nf)
  for (i in seq_len(nr)) {
    Ri <- ref[[i]]
    for (j in seq_len(nf)) {
      Fj <- found[[j]]
      ov <- length(intersect(Ri, Fj))
      if (ov == 0L) next
      p <- ov / length(Fj); r <- ov / length(Ri)
      P[i, j] <- p; R[i, j] <- r
      F1[i, j] <- 2 * p * r / (p + r)
      J[i, j] <- ov / length(union(Ri, Fj))
    }
  }
  list(f1 = F1, precision = P, recall = R, jaccard = J)
}

# Pad a rectangular reward matrix to square with zero (dummy) rows/cols.
.pad_square <- function(M) {
  nr <- nrow(M); nf <- ncol(M)
  n <- max(nr, nf)
  if (nr == nf) return(list(M = M, nr = nr, nf = nf, n = n))
  Mp <- matrix(0, n, n)
  Mp[seq_len(nr), seq_len(nf)] <- M
  list(M = Mp, nr = nr, nf = nf, n = n)
}

#' Align found modules to reference regulons
#'
#' @description
#' \code{$new(reference)} stores reference modules (e.g. RegulonDB regulons);
#' \code{$set_found(modules)} stores the modules under test;
#' \code{$align()} aligns each reference module to the found modules and
#' returns one row per reference module with expected precision/recall/F1;
#' \code{$summary()} aggregates into headline numbers. The active bindings
#' expose the assignment state: \code{$alignment} / \code{$assignment} are
#' the stable Hungarian headline, \code{$soft_alignment} /
#' \code{$soft_assignment} the soft uncertainty layer, and
#' \code{$soft_uncertainty} per-regulon dummy mass and entropy.
#' @export
ModuleTester <- R6::R6Class("ModuleTester",
  public = list(
    #' @description Store and validate reference modules.
    #' @param reference named list of character vectors (gene IDs). Names are
    #'   regulon/module labels; duplicates within a module are dropped.
    initialize = function(reference = NULL) {
      if (!is.null(reference)) self$set_reference(reference)
      invisible(self)
    },

    #' @description Replace the reference module set.
    #' @param reference named list of character vectors.
    set_reference = function(reference) {
      private$reference_ <- private$.validate_modules(reference, "reference")
      private$cache_ <- list()
      private$soft_epsilon_ <- 0.1
      invisible(self)
    },

    #' @description Store the found modules to be scored.
    #' @param modules list of character vectors (names optional).
    set_found = function(modules) {
      private$found_ <- private$.validate_modules(modules, "found",
                                                  named = FALSE)
      private$cache_ <- list()
      private$soft_epsilon_ <- 0.1
      invisible(self)
    },

    #' @description Align reference modules to found modules.
    #'
    #' @param method "hungarian" (default): hard one-to-one assignment
    #'   maximizing total F1. This is the headline recovery metric: an
    #'   exactly recovered module set scores exactly 1, and one-to-one
    #'   matching prevents a single found "hairball" from earning credit for
    #'   several regulons. "soft": entropic-regularized probabilistic
    #'   assignment as an uncertainty/ambiguity layer; reported scores are
    #'   expectations under each regulon's assignment distribution.
    #'   "best": independent argmax per regulon (many-to-one allowed).
    #' @param epsilon temperature for "soft": smaller concentrates mass on
    #'   the best matches (limit: Hungarian), larger spreads it. Soft
    #'   results are cached per epsilon value, so re-running with a new
    #'   temperature never silently reuses a stale one.
    #' @return data.frame, one row per reference module: reference, size_ref,
    #'   f1, precision, recall, jaccard (expectations under "soft", exact
    #'   under the hard methods), top_found (highest-probability / assigned
    #'   found module), top_prob (its probability; 1 under hard methods).
    align = function(method = c("hungarian", "soft", "best"),
                     epsilon = 0.1) {
      method <- match.arg(method)
      if (method == "soft") {
        if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0)
          stop("epsilon must be a single finite positive number",
               call. = FALSE)
        private$soft_epsilon_ <- epsilon
      }
      key <- private$.cache_key(method, epsilon)
      res <- private$.compute_alignment(method, epsilon)
      # Per-method (and per-epsilon) cache: headline ("hungarian") and
      # uncertainty ("soft") state coexist, so requesting one never
      # overwrites the other.
      private$cache_[[key]] <- res
      res$alignment
    },

    #' @description Aggregate the alignment into headline numbers.
    #'
    #' Headline aggregates default to "hungarian" assignment (exact recovery
    #' scores exactly 1); pass method = "soft" to aggregate the soft
    #' uncertainty layer instead. Each method's alignment is computed and
    #' cached independently (soft: per epsilon value), so a headline summary
    #' never overwrites a cached soft alignment.
    #'
    #' Returns a list: n_ref, n_found, method, epsilon (NA for hard
    #' methods), mean_f1 (unweighted mean of per-regulon F1), weighted_f1
    #' (regulon-size-weighted mean), median_f1, recovery_rate (fraction of
    #' regulons with F1 >= 0.5), mean precision / recall; for
    #' method = "soft" additionally mean_dummy_mass, mean_entropy, converged
    #' and sinkhorn_iters from the cached soft fit.
    summary = function(method = c("hungarian", "soft", "best"),
                       epsilon = 0.1) {
      method <- match.arg(method)
      key <- private$.cache_key(method, epsilon)
      if (is.null(private$cache_[[key]]))
        private$cache_[[key]] <- private$.compute_alignment(method,
                                                            epsilon)
      e <- private$cache_[[key]]
      al <- e$alignment
      w <- al$size_ref / sum(al$size_ref)
      out <- list(
        n_ref = nrow(al),
        n_found = length(private$found_),
        method = method,
        epsilon = e$epsilon,
        mean_f1 = mean(al$f1),
        weighted_f1 = sum(w * al$f1),
        median_f1 = stats::median(al$f1),
        recovery_rate = mean(al$f1 >= 0.5),
        mean_precision = mean(al$precision),
        mean_recall = mean(al$recall)
      )
      if (method == "soft") {
        out$mean_dummy_mass <- mean(e$dummy_mass)
        out$mean_entropy <- mean(e$entropy)
        out$converged <- e$converged
        out$sinkhorn_iters <- e$iters
      }
      out
    },

    #' @description Summarize tester state.
    print = function(...) {
      nr <- if (is.null(private$reference_)) 0L else length(private$reference_)
      nf <- if (is.null(private$found_)) 0L else length(private$found_)
      cat("<ModuleTester>", nr, "reference modules,", nf, "found modules")
      # Headline only: never trigger a fresh computation from print(), and
      # never let soft state leak into the headline display.
      if (!is.null(private$cache_[["hungarian"]])) {
        s <- self$summary()
        cat(sprintf(" [%s], mean F1 %.3f", s$method, s$mean_f1))
      }
      cat("\n")
      invisible(self)
    }
  ),

  active = list(
    #' @field reference the validated reference module list.
    reference = function() private$.need(private$reference_,
                                         "set_reference()"),
    #' @field found the validated found-module list.
    found = function() private$.need(private$found_, "set_found()"),
    #' @field alignment the Hungarian headline alignment table (computes
    #'   and caches it if needed). This is the stable headline state: it
    #'   never reflects a soft or best alignment.
    alignment = function() private$.cached("hungarian", 0.1)$alignment,
    #' @field assignment the Hungarian headline assignment matrix,
    #'   reference x found, 0/1 (computes and caches if needed).
    assignment = function() private$.cached("hungarian", 0.1)$P,
    #' @field best_alignment the independent-best alignment table
    #'   (many-to-one allowed; diagnostic, computes and caches if needed).
    best_alignment = function() private$.cached("best", 0.1)$alignment,
    #' @field best_assignment the independent-best assignment matrix,
    #'   reference x found, 0/1 (computes and caches if needed).
    best_assignment = function() private$.cached("best", 0.1)$P,
    #' @field soft_alignment the soft alignment table at the current soft
    #'   epsilon (see $align(method = "soft")); expectations under each
    #'   regulon's assignment distribution. Computes and caches if needed.
    soft_alignment = function() {
      private$.cached("soft", private$soft_epsilon_)$alignment
    },
    #' @field soft_assignment the soft assignment matrix, reference x
    #'   found, at the current soft epsilon: row i is regulon i's
    #'   distribution over real found modules (rows sum to <= 1; the
    #'   shortfall is dummy/"unexplained" mass).
    soft_assignment = function() {
      private$.cached("soft", private$soft_epsilon_)$P
    },
    #' @field soft_uncertainty per-regulon uncertainty at the current soft
    #'   epsilon: data.frame with reference, size_ref, dummy_mass
    #'   (probability assigned to dummy "unexplained" modules), entropy
    #'   (nats; over the full row distribution *including* dummy mass, so
    #'   unexplained regulons count as uncertain), top_found and top_prob.
    soft_uncertainty = function() {
      e <- private$.cached("soft", private$soft_epsilon_)
      data.frame(reference = e$alignment$reference,
                 size_ref = e$alignment$size_ref,
                 dummy_mass = e$dummy_mass,
                 entropy = e$entropy,
                 top_found = e$alignment$top_found,
                 top_prob = e$alignment$top_prob,
                 stringsAsFactors = FALSE)
    }
  ),

  private = list(
    reference_ = NULL, found_ = NULL, cache_ = list(),
    # Epsilon used by the soft_* active bindings: the most recent epsilon
    # passed to $align(method = "soft"), or the default 0.1.
    soft_epsilon_ = 0.1,

    # Cache key for one alignment: hard methods by name, soft per epsilon
    # so different temperatures never alias each other.
    .cache_key = function(method, epsilon) {
      if (method == "soft") paste0("soft:", as.character(epsilon))
      else method
    },

    # Fetch (computing and caching on miss) one alignment entry.
    .cached = function(method, epsilon) {
      key <- private$.cache_key(method, epsilon)
      if (is.null(private$cache_[[key]]))
        private$cache_[[key]] <- private$.compute_alignment(method,
                                                            epsilon)
      private$cache_[[key]]
    },

    # Compute (but do not cache) one alignment; returns a list with
    # alignment (data.frame), P (nr x nf assignment matrix), method,
    # epsilon (NA for hard methods), dummy_mass and entropy (soft only;
    # NA otherwise), and convergence diagnostics (soft only).
    .compute_alignment = function(method, epsilon) {
      ref <- private$.need(private$reference_, "set_reference()")
      found <- private$.need(private$found_, "set_found()")
      st <- .pairwise_set_stats(ref, found)
      nr <- length(ref); nf <- length(found)
      out <- list(method = method,
                  epsilon = if (method == "soft") epsilon else NA_real_,
                  dummy_mass = rep(NA_real_, nr),
                  entropy = rep(NA_real_, nr),
                  converged = NA,
                  iters = NA_integer_,
                  marginal_error = NA_real_)

      if (method == "best") {
        # Independent argmax per reference module (many-to-one allowed).
        P <- matrix(0, nr, nf)
        for (i in seq_len(nr)) {
          j <- which.max(st$f1[i, ])
          if (st$f1[i, j] > 0) P[i, j] <- 1
        }
        rows <- lapply(seq_len(nr), function(i) {
          j <- which(P[i, ] > 0)
          private$.row(names(ref)[i], length(ref[[i]]), st, i, j,
                       names(found), prob = 1)
        })
      } else {
        pad <- .pad_square(st$f1)
        n <- pad$n
        if (method == "hungarian") {
          # Maximize total F1 <=> minimize -F1.
          a <- .hungarian(-pad$M)   # a[col] = row
          Ppad <- matrix(0, n, n)
          for (j in seq_len(n)) Ppad[a[j], j] <- 1
        } else {
          if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon <= 0)
            stop("epsilon must be a single finite positive number",
                 call. = FALSE)
          sk <- .sinkhorn(pad$M, epsilon = epsilon)
          Ppad <- sk$P
          out$converged <- sk$converged
          out$iters <- sk$iters
          out$marginal_error <- sk$marginal_error
        }
        # Real x real submatrix, row-normalized: each row is the regulon's
        # assignment distribution over real found modules (rows sum to <= 1;
        # the shortfall is mass on dummy modules, i.e. "unexplained", and it
        # penalizes the expectation below).
        P <- Ppad[seq_len(nr), seq_len(nf), drop = FALSE]
        # Invariant: Ppad has uniform 1/n row margins, so rowmass == 1/n > 0
        # for every real row (soft) or == 1 (hungarian permutation rows).
        rowmass <- rowSums(Ppad[seq_len(nr), , drop = FALSE])
        P <- P / rowmass
        rows <- lapply(seq_len(nr), function(i) {
          private$.expected_row(names(ref)[i], length(ref[[i]]), st, i,
                                P[i, ], names(found))
        })
      }
      if (method == "soft") {
        # Dummy mass: probability the regulon is explained by no found
        # module. Entropy (nats) is over the full row distribution
        # *including* dummy mass, so unexplained regulons count as
        # uncertain rather than as confident misses.
        dm <- 1 - rowSums(P)
        dm <- pmax(dm, 0)
        out$dummy_mass <- dm
        out$entropy <- vapply(seq_len(nr), function(i) {
          q <- c(P[i, ], dm[i])
          q <- q[q > 0]
          -sum(q * log(q))
        }, numeric(1))
      }
      al <- do.call(rbind, rows)
      rownames(al) <- NULL
      out$alignment <- al
      out$P <- P
      out
    },

    .need = function(value, stage) {
      if (is.null(value))
        stop("not available yet: run $", stage, " first", call. = FALSE)
      value
    },

    # One alignment-table row for a hard assignment of regulon i to found j.
    .row = function(rn, size_ref, st, i, j, found_names, prob) {
      if (!length(j)) {
        return(data.frame(reference = rn, size_ref = size_ref,
                          f1 = 0, precision = 0, recall = 0, jaccard = 0,
                          top_found = NA_character_, top_prob = 0,
                          stringsAsFactors = FALSE))
      }
      data.frame(reference = rn, size_ref = size_ref,
                 f1 = st$f1[i, j], precision = st$precision[i, j],
                 recall = st$recall[i, j], jaccard = st$jaccard[i, j],
                 top_found = found_names[j], top_prob = prob,
                 stringsAsFactors = FALSE)
    },

    # One alignment-table row: expectations under assignment weights w.
    # w is already a distribution over real found modules (sums to <= 1);
    # the shortfall is mass on dummy ("unexplained") modules and it stays
    # out of the expectation, i.e. it penalizes the score.
    .expected_row = function(rn, size_ref, st, i, w, found_names) {
      if (sum(w) <= 0) {
        return(data.frame(reference = rn, size_ref = size_ref,
                           f1 = 0, precision = 0, recall = 0, jaccard = 0,
                           top_found = NA_character_, top_prob = 0,
                           stringsAsFactors = FALSE))
      }
      j <- which.max(w)
      data.frame(reference = rn, size_ref = size_ref,
                 f1 = sum(w * st$f1[i, ]),
                 precision = sum(w * st$precision[i, ]),
                 recall = sum(w * st$recall[i, ]),
                 jaccard = sum(w * st$jaccard[i, ]),
                 top_found = found_names[j], top_prob = w[j],
                 stringsAsFactors = FALSE)
    },

    # Coerce to a named list of unique non-empty character vectors.
    .validate_modules = function(modules, what, named = TRUE) {
      if (!is.list(modules) || length(modules) == 0L)
        stop(what, " modules must be a non-empty list of gene-ID vectors")
      nm <- names(modules)
      if (named && (is.null(nm) || any(!nzchar(nm))))
        stop(what, " modules must be a *named* list (regulon labels)")
      if (is.null(nm)) nm <- paste0("module_", seq_along(modules))
      out <- lapply(seq_along(modules), function(i) {
        v <- modules[[i]]
        if (is.factor(v)) v <- as.character(v)
        if (!is.character(v) && !is.numeric(v) && !is.integer(v))
          stop(what, " module '", nm[i], "' must be a vector of gene IDs")
        v <- unique(as.character(v))
        v <- v[!is.na(v) & nzchar(v)]
        if (!length(v)) stop(what, " module '", nm[i], "' is empty")
        v
      })
      names(out) <- nm
      out
    }
  )
)

#' Find modules as connected components of the attention graph
#'
#' Symmetrizes a thresholded attention operator (\code{(A + t(A)) / 2 > 0})
#' and returns its connected components as gene-ID lists. This is the
#' toy-grade module finder: on real data the components of a sparse attention
#' graph are a reasonable first cut, but richer community detection is a
#' separate research question. Isolated genes are returned as singletons so
#' that no gene silently disappears from the evaluation.
#'
#' @param operator square numeric matrix (genes x genes), e.g.
#'   \code{ClrAttention$operator}; dimnames used as gene IDs when present.
#' @param min_size drop components smaller than this (default 1 keeps all).
#' @return named list of character vectors, largest component first.
#' @export
find_modules <- function(operator, min_size = 1L) {
  if (!is.matrix(operator) || !is.numeric(operator) ||
      nrow(operator) != ncol(operator))
    stop("operator must be a square numeric matrix")
  G <- nrow(operator)
  ids <- rownames(operator)
  if (is.null(ids)) ids <- as.character(seq_len(G))
  adj <- ((operator + t(operator)) / 2) > 0
  diag(adj) <- FALSE
  # BFS connected components.
  comp <- integer(G)
  cid <- 0L
  for (s in seq_len(G)) {
    if (comp[s] != 0L) next
    cid <- cid + 1L
    queue <- s
    comp[s] <- cid
    while (length(queue)) {
      v <- queue[1L]
      queue <- queue[-1L]
      nbrs <- which(adj[v, ] & comp == 0L)
      comp[nbrs] <- cid
      queue <- c(queue, nbrs)
    }
  }
  mods <- split(ids, comp)
  mods <- mods[order(lengths(mods), decreasing = TRUE)]
  names(mods) <- paste0("module_", seq_along(mods))
  mods <- mods[lengths(mods) >= min_size]
  if (!length(mods)) stop("no components of size >= min_size")
  mods
}
