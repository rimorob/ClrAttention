# EXPLORATORY (docs/plans/explore_stack_plan.md). Stacked a3r layers: generic layer function, readouts, controls, label-free diagnostics.
# Sourced by stack_k562.R and stack_m3d.R after robust_z_lib.R (rz_operator, rz_lazy_P, rz_stationary). Not for the paper or the pre-read.

# ---- readouts ------------------------------------------------------------------------------------------------
# a3r: bilateral CLR on log P^t, row statistics over reachable pairs only (docs/plans/corum_analysis.md reference code; identical to
# ideas47_lib.R::a3r). Pt: dense P^t. Returns the symmetric z-scale matrix, truncated at 0, zero diagonal.
sk_a3r <- function(Pt) {
  L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  Z <- (L - mu) / sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1)); Z[is.na(Z) | Z < 0] <- 0
  S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S
}
# c2_z: a3r-style bilateral row standardisation applied again to S directly (no log, no softmax, no walk): row mean and SD over the
# off-diagonal entries of S (zeros included), negatives and non-finite z set to 0, Stouffer combination, zero diagonal.
sk_c2z <- function(S) {
  diag(S) <- NA; mu <- rowMeans(S, na.rm = TRUE); n <- rowSums(!is.na(S))
  Z <- (S - mu) / sqrt(rowSums((S - mu)^2, na.rm = TRUE) / (n - 1)); Z[!is.finite(Z) | Z < 0] <- 0
  Z <- (Z + t(Z)) / sqrt(2); diag(Z) <- 0; Z
}
# c2_log: the literal a3r function applied to S (log S, zeros unreachable)
sk_c2log <- function(S) { S[S <= 0] <- 0; sk_a3r(S) }

# reachability of a dense P^t: share of off-diagonal pairs with P^t = 0, rows with fewer than 2 reachable off-diagonal entries
sk_reach <- function(Pt) {
  L <- Pt; diag(L) <- NA; ok <- is.finite(L) & L > 0 & !is.na(L); G <- nrow(L)
  c(unreach_share = 1 - sum(ok) / (G * (G - 1)), rows_lt2_reachable = sum(rowSums(ok) < 2L))
}

# ---- one stacked layer ---------------------------------------------------------------------------------------
# Sin: symmetric score matrix (CLR for layer 1, S_k for layer k + 1); or A: a prebuilt operator (reused for the matched-depth control).
# ts_out: depths at which the a3r readout is returned; t_diag: depth at which reachability is recorded.
# Operator = package path (rz_operator: top-50 positive candidates per row, keff-10 softmax, row normalisation, self-loop for empty rows),
# lazy walk 0.5, P^t by repeated multiplication (P sparse).
sk_layer <- function(Sin = NULL, A = NULL, ts_out = 9L, t_diag = 9L, alpha = 0.5) {
  if (is.null(A)) A <- rz_operator(Sin)
  P <- rz_lazy_P(A, alpha); outs <- list(); reach <- list(); Pt <- NULL
  for (t in seq_len(max(ts_out))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (t %in% ts_out) { outs[[as.character(t)]] <- sk_a3r(Pt); reach[[as.character(t)]] <- sk_reach(Pt) }
    if (t %% 8L == 0L) invisible(gc(verbose = FALSE))
  }
  rm(Pt); invisible(gc(verbose = FALSE))
  list(A = A, P = P, outs = outs, reach = reach)
}

# ---- label-free diagnostics per layer ------------------------------------------------------------------------
# k: layer index (1 = a3r of the CLR operator); Sin: the score the layer's operator was built from; A, P: the layer operator and walk;
# reach: sk_reach at the diagnostic depth; Sout: the layer's a3r output; top1_prev / top1_clr: logical top-1% pair sets (upper triangle) of
# S_{k-1} and of the CLR score; mi_rm: MI row mean; pi1: stationary distribution of layer 1; UT: upper-triangle indices.
sk_top1 <- function(S, UT) { v <- S[UT]; q <- unname(stats::quantile(v, 0.99)); v >= q & v > 0 }
sk_jaccard <- function(a, b) sum(a & b) / sum(a | b)
sk_diag <- function(dataset, k, Sin, A, P, reach, Sout, top1_prev, top1_clr, mi_rm, pi1, UT, t = 9L) {
  G <- nrow(A); selfloop <- diag(A) > 0
  cnt <- rowSums(A > 0) - selfloop; npos <- rowSums(Sin > 0)
  pik <- rz_stationary(P); cs <- colSums(A); top <- order(cs, decreasing = TRUE)[seq_len(max(1L, round(0.01 * G)))]
  t1 <- sk_top1(Sout, UT)
  data.frame(dataset = dataset, layer = k, t = t,
    rows_no_candidate = sum(npos == 0L), frac_rows_no_candidate = mean(npos == 0L),
    rows_lt10_candidates = sum(cnt < 10L), median_candidates = stats::median(cnt[!selfloop]),
    median_eff_support = stats::median(1 / rowSums(A[!selfloop, , drop = FALSE]^2)),
    unreach_share = reach[["unreach_share"]], rows_lt2_reachable = reach[["rows_lt2_reachable"]],
    spearman_pi_mi = stats::cor(pik, mi_rm, method = "spearman"), spearman_pi_pi1 = stats::cor(pik, pi1, method = "spearman"),
    pi_max_over_median = max(pik) / stats::median(pik), inweight_share_top1pct = sum(cs[top]) / sum(cs),
    S_frac_nonzero = mean(Sout[UT] > 0), S_rows_all_zero = sum(rowSums(Sout > 0) == 0L),
    spearman_Srowmean_mi = stats::cor(rowMeans(Sout), mi_rm, method = "spearman"),
    jaccard_top1_prev = sk_jaccard(t1, top1_prev), jaccard_top1_clr = sk_jaccard(t1, top1_clr), stringsAsFactors = FALSE)
}

# ---- scoring (used by the M3D driver; K562 uses ideas47_lib.R::score_readout) ---------------------------------
# bench: list(truth = list(pairs = comembership_pairs...)); halves: list(truth = list of 100 gene index vectors).
# Returns a function score(name, S) giving list(name, raw, hub_att, hub_mi), each a per-truth list(aupr, jack); cached as readout__<name>.rds.
sk_make_scorer <- function(bench, halves, G, H_att, H_mi, nc, out, say) {
  jack_one <- function(sc, b, keeps) { i <- b$pairs$i; j <- b$pairs$j; lab_ <- b$pairs$label
    unlist(parallel::mclapply(keeps, function(keep) { inh <- logical(G); inh[keep] <- TRUE; ix <- which(inh[i] & inh[j]); ap_fast(sc[ix], lab_[ix]) }, mc.cores = nc)) }
  pooled <- function(S) { invisible(gc(verbose = FALSE)); lapply(stats::setNames(names(bench), names(bench)), function(e) { b <- bench[[e]]
    sc <- S[cbind(b$pairs$i, b$pairs$j)]; list(aupr = ap_fast(sc, b$pairs$label), jack = jack_one(sc, b, halves[[e]])) }) }
  function(name, S) {
    cf <- file.path(out, "cache", paste0("readout__", name, ".rds"))
    if (file.exists(cf)) { say("cached: ", name); return(invisible(readRDS(cf))) }
    t1 <- Sys.time(); S[!is.finite(S)] <- 0; diag(S) <- 0
    res <- list(name = name, raw = pooled(S)); invisible(gc(verbose = FALSE))
    Sh <- rz_hub_stratify(S, H_att); res$hub_att <- pooled(Sh); rm(Sh); invisible(gc(verbose = FALSE))
    Sh <- rz_hub_stratify(S, H_mi); res$hub_mi <- pooled(Sh); rm(Sh); invisible(gc(verbose = FALSE))
    saveRDS(res, cf); say(sprintf("scored %-12s (%.0f s)", name, as.numeric(difftime(Sys.time(), t1, units = "secs")))); invisible(res)
  }
}

# stored-gate bookkeeping shared by both drivers: rows of (check, max_abs_diff, pass)
sk_gate_add <- function(rows, check, value, tol = 1e-9) { rows[[length(rows) + 1L]] <- data.frame(check = check, max_abs_diff = value, pass = is.finite(value) && value < tol); rows }
sk_gate_finish <- function(rows, file, say) {
  g <- do.call(rbind, rows); utils::write.csv(g, file, row.names = FALSE); print(g, row.names = FALSE)
  if (!all(g$pass)) stop("EQUIVALENCE GATE FAILED (", file, "): no stacked readout is computed or reported")
  say("EQUIVALENCE GATE PASSED: ", basename(file)); invisible(g)
}
