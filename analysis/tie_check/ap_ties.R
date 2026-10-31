# Average precision (pooled AUPR) under different tie conventions.
#
# analysis/regulondb.R::ap_fast gives every pair in a block of exactly tied
# scores the precision at the END of the block (a positive inside a tie is
# credited as if all of the block's negatives ranked above it). This file
# computes, from the same sorted blocks, four conventions:
#   old         = ap_fast's convention (reproduced exactly)
#   expected    = expected AP when each tied block is ordered uniformly at random
#                 (exact analytic expectation, no Monte Carlo)
#   optimistic  = positives first within every block (upper bound)
#   negfirst    = negatives first within every block (strict lower bound; the
#                 j-th positive of a block gets the precision at its own rank)
#
# For a block with m pairs, p positives, preceded by k0 pairs holding tp0
# positives, AP numerator contributions (sum of precision at its positives):
#   old        p (tp0 + p) / (k0 + m)
#   optimistic sum_{j=1..p} (tp0 + j) / (k0 + j)
#   negfirst   sum_{j=1..p} (tp0 + j) / (k0 + m - p + j)
#   expected   (p/m) sum_{r=1..m} (tp0 + 1 + (r-1)(p-1)/(m-1)) / (k0 + r)
# The expected form follows from linearity: a positive at block position r has
# a Hypergeometric(m-1, p-1, r-1) number of earlier positives in the block,
# mean (r-1)(p-1)/(m-1), and precision is linear in that count. All sums are
# of the form H(k, n) = sum_{r=1..n} 1/(k + r), evaluated by digamma
# differences for long blocks and by direct summation for short ones.
# Missing scores are treated as -Inf (one tied block at the bottom), as in ap_fast.

H_sum <- function(k, n) {          # sum_{r=1..n} 1 / (k + r), vectorised; n = 0 gives 0
  out <- numeric(length(k)); big <- n > 32
  if (any(big)) out[big] <- digamma(k[big] + n[big] + 1) - digamma(k[big] + 1)
  sm <- which(!big & n > 0)
  if (length(sm)) {
    idx <- rep.int(sm, n[sm]); r <- sequence(n[sm])
    out[sm] <- as.numeric(rowsum(1 / (k[idx] + r), match(idx, sm), reorder = TRUE))
  }
  out
}

# AP numerators from block summaries (vectors over blocks), divided by n1
ap_from_blocks <- function(k0, m, tp0, p, n1) {
  k0 <- as.numeric(k0); m <- as.numeric(m); tp0 <- as.numeric(tp0); p <- as.numeric(p)
  hp <- p > 0
  k0 <- k0[hp]; m <- m[hp]; tp0 <- tp0[hp]; p <- p[hp]
  b <- ifelse(m > 1, (p - 1) / pmax(m - 1, 1), 0)
  kn <- k0 + m - p
  c(old        = sum(p * (tp0 + p) / (k0 + m)) / n1,
    expected   = sum(p / m * (b * m + (tp0 + 1 - b * (k0 + 1)) * H_sum(k0, m))) / n1,
    optimistic = sum(p + (tp0 - k0) * H_sum(k0, p)) / n1,
    negfirst   = sum(p + (tp0 - kn) * H_sum(kn, p)) / n1)
}

# Block summaries of a (sorted) label vector with block end markers
block_stats <- function(l, last) {
  k <- which(last); cl <- cumsum(l)[last]
  m <- diff(c(0, k)); p <- diff(c(0, cl))
  list(k0 = k - m, m = m, tp0 = cl - p, p = p)
}

# All four conventions for one score/label vector.
ap_ties <- function(score, label) {
  score[is.na(score)] <- -Inf
  o <- order(score, decreasing = TRUE, method = "radix")
  s <- score[o]; l <- label[o]
  last <- c(s[-1L] != s[-length(s)], TRUE)
  b <- block_stats(l, last)
  ap_from_blocks(b$k0, b$m, b$tp0, b$p, sum(label))
}

# Tie summary: fraction of pairs / positives inside blocks of size > 1, size of
# the largest block, number of distinct values
tie_summary <- function(score, label) {
  score[is.na(score)] <- -Inf
  o <- order(score, decreasing = TRUE, method = "radix")
  s <- score[o]; l <- label[o]
  last <- c(s[-1L] != s[-length(s)], TRUE)
  b <- block_stats(l, last); t <- b$m > 1
  c(n_pairs = length(score), n_pos = sum(label), n_blocks = length(b$m),
    n_tied_pairs = sum(b$m[t]), n_tied_pos = sum(b$p[t]), max_block = max(b$m),
    max_block_pos = b$p[which.max(b$m)])
}

# Prepared version for the delete-half jackknife: sort once, then restrict the
# sorted order to a subset of pairs (a subset of a sorted vector stays sorted,
# and tied blocks stay contiguous), so each replicate costs one pass, no sort.
ap_prepare <- function(score, label) {
  score[is.na(score)] <- -Inf
  o <- order(score, decreasing = TRUE, method = "radix")
  s <- score[o]
  list(o = o, g = cumsum(c(TRUE, s[-1L] != s[-length(s)])), l = label[o])
}
ap_ties_subset <- function(pp, ix_mask) {
  sel <- ix_mask[pp$o]
  g <- pp$g[sel]; l <- pp$l[sel]
  last <- c(g[-1L] != g[-length(g)], TRUE)
  b <- block_stats(l, last)
  ap_from_blocks(b$k0, b$m, b$tp0, b$p, sum(l))
}
