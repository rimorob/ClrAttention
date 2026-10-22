# Gene filters for building a CLR context on single-cell data (sourced by
# beeline.R, beeline_boot.R and run directly for counts).
#
# Ranking genes by variance is a poor context filter on scRNA-seq: it favours
# noisy, low-expressed genes and drops well-expressed genes whose variation is
# real. Instead keep genes that are DRIVEN by something in the cell:
#   multimodal  Hartigan's dip test on log expression (p < dip_p); sufficient
#               but not necessary;
#   pc-driven   the gene's expression depends on cell state: in a forward
#               stepwise regression of the gene on the top K principal
#               components of the cells (K by Horn's parallel analysis), at
#               least one PC enters. Because PCs are orthogonal, the first
#               forward step of stepAIC has a closed form: some PC enters iff
#               n * log(1 / (1 - r_max^2)) > penalty, with r_max the gene's
#               largest |correlation| with a PC and penalty 2 (AIC) or log(n)
#               (BIC). (Checked against MASS::stepAIC on 60 genes: 60/60 agree.)
# Both criteria detect dependence on SOMETHING, and an uncorrected batch is
# such a something: run them only on batch-corrected data. BEELINE datasets
# are single studies without batch labels, so they are used as they are.

# Horn's parallel analysis on genes x cells: number of components whose
# eigenvalue exceeds the q-quantile of eigenvalues from data with every gene
# permuted independently across cells. Returns K and the cell scores.
horn_pcs <- function(X, reps = 5L, q = 0.95, kmax = 50L) {
  Z <- X - rowMeans(X); s <- sqrt(rowSums(Z^2) / (ncol(Z) - 1)); s[s <= 0] <- 1; Z <- Z / s
  ev <- function(M) { e <- eigen(crossprod(M), symmetric = TRUE, only.values = TRUE)$values
                      e[seq_len(kmax)] / (nrow(M) - 1) }
  obs <- ev(Z)
  null <- vapply(seq_len(reps), function(r) ev(t(apply(Z, 1, sample))), numeric(kmax))
  thr <- apply(null, 1, stats::quantile, q)
  K <- max(1L, which(obs <= thr)[1] - 1L, na.rm = TRUE)
  if (is.na(K) || K < 1L) K <- kmax
  e <- eigen(crossprod(Z), symmetric = TRUE)
  list(K = K, scores = e$vectors[, seq_len(K), drop = FALSE], eig = obs, null_thr = thr)
}

# Per-gene filter statistics for an expression matrix (genes x cells, log scale).
gene_filter_stats <- function(X, reps = 5L, dip_p = 0.01, kmax = 50L) {
  hp <- horn_pcs(X, reps = reps, kmax = kmax)
  n <- ncol(X)
  r <- abs(stats::cor(t(X), hp$scores)); r[is.na(r)] <- 0
  r2 <- apply(r, 1, max)^2
  gain <- n * log(1 / pmax(1 - r2, 1e-12))
  dip <- vapply(seq_len(nrow(X)), function(g) diptest::dip.test(X[g, ])$p.value, 0)
  data.frame(gene = rownames(X), mean = rowMeans(X), frac_nonzero = rowMeans(X > 0),
             variance = apply(X, 1, stats::var), r2_max = r2,
             pc_aic = gain > 2, pc_bic = gain > log(n), dip_p = dip, multimodal = dip < dip_p,
             K = hp$K, stringsAsFactors = FALSE)
}

# Context genes: expressed genes that are multimodal or PC-driven (AIC or BIC).
driven_genes <- function(X, rule = c("aic", "bic"), ...) {
  rule <- match.arg(rule)
  st <- gene_filter_stats(X, ...)
  keep <- st$multimodal | (if (rule == "aic") st$pc_aic else st$pc_bic)
  st$gene[keep]
}
