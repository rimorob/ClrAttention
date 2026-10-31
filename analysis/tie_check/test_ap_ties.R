#!/usr/bin/env Rscript
# Tests for analysis/tie_check/ap_ties.R. Run from the repo root.
source("analysis/regulondb.R"); source("analysis/tie_check/ap_ties.R")
ok <- function(cond, msg) { cat(if (isTRUE(cond)) "PASS " else "FAIL ", msg, "\n"); if (!isTRUE(cond)) quit(status = 1) }

# brute force over every arrangement of positives within tied blocks
brute <- function(score, label) {
  blocks <- split(seq_along(score), factor(score, levels = sort(unique(score), decreasing = TRUE)))
  n1 <- sum(label); N <- length(score)
  ap_of <- function(lab_sorted) { tp <- cumsum(lab_sorted); sum((tp / seq_along(tp))[lab_sorted == 1]) / n1 }
  # enumerate all distinct within-block orderings of labels
  perms <- function(v) { if (length(v) <= 1) return(list(v)); u <- unique(v)
    do.call(c, lapply(u, function(x) { rest <- v; rest <- rest[-match(x, rest)]; lapply(perms(rest), function(r) c(x, r)) })) }
  opts <- lapply(blocks, function(ix) perms(label[ix]))
  grid <- expand.grid(lapply(opts, seq_along)); vals <- numeric(nrow(grid))
  for (g in seq_len(nrow(grid))) vals[g] <- ap_of(unlist(mapply(function(o, i) o[[i]], opts, as.list(grid[g, ]), SIMPLIFY = FALSE)))
  # arrangements are enumerated by distinct label sequences; weight by multiplicity (uniform over permutations of items)
  w <- vapply(seq_len(nrow(grid)), function(g) prod(mapply(function(o, i, ix) {
    lab <- o[[i]]; factorial(sum(lab == 1)) * factorial(sum(lab == 0)) / 1 }, opts, as.list(grid[g, ]), blocks)), 0)
  ap_ng <- vapply(seq_len(nrow(grid)), function(g) ap_of(unlist(mapply(function(o, i) o[[i]], opts, as.list(grid[g, ]), SIMPLIFY = FALSE))), 0)
  c(expected = sum(w * ap_ng) / sum(w), optimistic = max(ap_ng), negfirst = min(ap_ng))
}

set.seed(1)
# 1. no ties: all conventions equal ap_fast
for (rep in 1:20) { n <- sample(5:400, 1); s <- rnorm(n); l <- rbinom(n, 1, runif(1, .05, .6)); if (!sum(l)) l[1] <- 1L
  a <- ap_ties(s, l); ok(all(abs(a - ap_fast(s, l)) < 1e-12), sprintf("no ties n=%d: all four = ap_fast", n)) }
# 2. old convention reproduces ap_fast with ties (incl. NA scores)
for (rep in 1:20) { n <- sample(5:2000, 1); s <- sample(c(0, 0, 0, 1:6, NA), n, TRUE) + 0; l <- rbinom(n, 1, runif(1, .05, .6)); if (!sum(l)) l[1] <- 1L
  a <- ap_ties(s, l); ok(abs(a[["old"]] - ap_fast(s, l)) < 1e-12, sprintf("ties n=%d: old = ap_fast", n))
  ok(a[["negfirst"]] <= a[["expected"]] + 1e-12 && a[["expected"]] <= a[["optimistic"]] + 1e-12, "  ordering negfirst <= expected <= optimistic") }
# 3. exact expectation against brute-force enumeration on small cases
cases <- list(c(1,1,1), c(1,2,2,2,3), c(1,1,2,2,3,3,3))
for (cs in cases) for (lab in list(c(1,0,1), c(0,1,0,1,1), c(1,0,0,1,0,1,1))[1:3]) {
  n <- length(cs); l <- rep_len(lab, n); if (sum(l) == 0) next
  b <- brute(cs, l); a <- ap_ties(cs, l)
  ok(all(abs(a[c("expected", "optimistic", "negfirst")] - b) < 1e-12), sprintf("brute force n=%d scores=%s labels=%s", n, paste(cs, collapse = ""), paste(l, collapse = "")))
}
for (rep in 1:30) { n <- sample(4:9, 1); s <- sample(1:3, n, TRUE); l <- rbinom(n, 1, .5); if (!sum(l)) l[1] <- 1L
  b <- brute(s, l); a <- ap_ties(s, l)
  ok(all(abs(a[c("expected", "optimistic", "negfirst")] - b) < 1e-12), sprintf("brute force random n=%d", n)) }
# 4. all tied: old = base rate, expected ~ base rate
l <- rbinom(1000, 1, .1); a <- ap_ties(rep(0, 1000), l)
ok(abs(a[["old"]] - mean(l)) < 1e-12, "all tied: old = base rate")
# 5. long blocks (digamma path) vs Monte Carlo with random tie-breaking
n <- 200000; s <- sample(c(0, 0, 0, 0, 1:50), n, TRUE); l <- rbinom(n, 1, ifelse(s > 40, .3, .05))
a <- ap_ties(s, l); mc <- replicate(200, { j <- runif(n); ap_fast(s + 0, l) * 0 + { o <- order(-s, j); t <- cumsum(l[o]); sum((t / seq_along(t))[l[o] == 1]) / sum(l) } })
ok(abs(a[["expected"]] - mean(mc)) < 4 * sd(mc) / sqrt(200) + 1e-6, sprintf("long blocks: analytic %.6f vs Monte Carlo %.6f +/- %.6f", a[["expected"]], mean(mc), sd(mc) / sqrt(200)))
# 6. subset path equals direct path
pp <- ap_prepare(s, l); m <- runif(n) < .5
ok(all(abs(ap_ties_subset(pp, m) - ap_ties(s[m], l[m])) < 1e-12), "prepared subset = direct")
# 7. real score vectors: tie_check.R checks the old-convention pooled and jackknife AUPR of every method/truth
#    against metrics.csv and jack_replicates.csv (summarize.R writes check_vs_existing_metrics.csv)
cat("all tests passed\n")
