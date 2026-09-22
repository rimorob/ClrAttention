test_that("exact module recovery gives headline F1 = 1", {
  ref <- list(regA = c("g1", "g2", "g3"), regB = c("g4", "g5"))
  mt <- ModuleTester$new(ref)
  mt$set_found(list(m1 = c("g1", "g2", "g3"), m2 = c("g4", "g5")))

  # Hungarian headline: exact recovery scores exactly 1.
  al <- mt$align(method = "hungarian")
  expect_equal(al$f1, c(1, 1))
  expect_equal(al$precision, c(1, 1))
  expect_equal(al$recall, c(1, 1))
  s <- mt$summary()
  expect_equal(s$method, "hungarian")
  expect_equal(s$mean_f1, 1)
  expect_equal(s$weighted_f1, 1)
  expect_equal(s$recovery_rate, 1)

  # Explicit headline bindings agree with align()/summary().
  expect_equal(mt$alignment$f1, c(1, 1))
  expect_equal(unname(mt$assignment), diag(2))
  expect_equal(mt$best_alignment$f1, c(1, 1))

  # Soft layer: entropy spreads a little mass, so expected F1 is just
  # below 1 -- but with no dummies there is no dummy mass and entropy is
  # tiny, not zero.
  mt$align(method = "soft", epsilon = 0.05)
  sal <- mt$soft_alignment
  expect_true(all(sal$f1 < 1))
  expect_true(all(sal$f1 > 0.99))
  u <- mt$soft_uncertainty
  expect_equal(u$dummy_mass, c(0, 0), tolerance = 1e-6)
  expect_true(all(u$entropy > 0))
  expect_true(all(u$entropy < 1e-3))

  # Headline state is untouched by the soft analysis.
  expect_equal(mt$alignment$f1, c(1, 1))
  expect_equal(mt$summary()$mean_f1, 1)
})

test_that("split modules lower recall, contaminated modules lower precision", {
  ref <- list(regA = c("g1", "g2", "g3", "g4"))
  mt <- ModuleTester$new(ref)
  # split: found modules each cover half -> recall 0.5, precision 1
  mt$set_found(list(m1 = c("g1", "g2"), m2 = c("g3", "g4")))
  al <- mt$align(method = "hungarian")
  expect_equal(al$recall, 0.5)
  expect_equal(al$precision, 1)
  expect_equal(al$f1, 2 * 1 * 0.5 / 1.5, tolerance = 1e-9)
  # contaminated: one found module covers all plus junk -> recall 1, precision 0.5
  mt$set_found(list(m1 = c("g1", "g2", "g3", "g4", "gx", "gy", "gz", "gw")))
  al <- mt$align(method = "hungarian")
  expect_equal(al$recall, 1)
  expect_equal(al$precision, 0.5)
  expect_equal(al$f1, 2 * 0.5 * 1 / 1.5, tolerance = 1e-9)
})

test_that("a missing regulon lands on dummy mass and scores 0", {
  ref <- list(regA = c("g1", "g2"), regB = c("g9", "g10"))
  mt <- ModuleTester$new(ref)
  mt$set_found(list(m1 = c("g1", "g2")))

  # Hard methods: regB is assigned to a dummy module -> F1 0.
  for (m in c("hungarian", "best")) {
    al <- mt$align(method = m)
    expect_equal(al$f1[al$reference == "regA"], 1, tolerance = 1e-9)
    expect_equal(al$f1[al$reference == "regB"], 0)
  }
  s <- mt$summary()
  expect_equal(s$mean_f1, 0.5)
  expect_equal(s$recovery_rate, 0.5)

  # Soft layer: regB's mass goes to the dummy ("unexplained"), so its
  # expected F1 is 0 and its dummy mass is ~1. Entropy includes the dummy
  # outcome, so a certainly-unexplained regulon has entropy 0.
  mt$align(method = "soft", epsilon = 0.1)
  sal <- mt$soft_alignment
  expect_equal(sal$f1[sal$reference == "regB"], 0, tolerance = 1e-9)
  expect_true(sal$f1[sal$reference == "regA"] > 0.99)
  u <- mt$soft_uncertainty
  expect_true(u$dummy_mass[u$reference == "regB"] > 0.99)
  expect_true(u$dummy_mass[u$reference == "regA"] < 0.01)
  # Entropy includes the dummy outcome, so a certainly-unexplained
  # regulon is nearly certain (tiny finite-temperature leakage only).
  expect_lt(u$entropy[u$reference == "regB"], 0.05)
})

test_that("one hairball cannot earn full credit for two regulons", {
  ref <- list(r1 = c("a", "b"), r2 = c("c", "d"))
  mt <- ModuleTester$new(ref)
  mt$set_found(list(hair = c("a", "b", "c", "d")))
  # F1 of either regulon against the hairball is 2/3; one-to-one matching
  # gives it to exactly one regulon.
  h <- mt$align(method = "hungarian")
  expect_equal(sum(h$f1 > 0), 1)
  expect_equal(sum(h$f1), 2 / 3, tolerance = 1e-9)
  # Independent-best (many-to-one) credits both -- the diagnostic contrast.
  b <- mt$align(method = "best")
  expect_equal(b$f1, c(2 / 3, 2 / 3), tolerance = 1e-9)
})

test_that("ambiguous matches split probability and raise entropy", {
  mt <- ModuleTester$new(list(r1 = c("a", "b", "c", "d")))
  mt$set_found(list(f1 = c("a", "b"), f2 = c("c", "d")))
  mt$align(method = "soft", epsilon = 0.5)
  P <- mt$soft_assignment
  # symmetric 2/3 vs 2/3 rewards -> mass splits evenly
  expect_equal(P[1, ], c(0.5, 0.5), tolerance = 1e-6)
  u <- mt$soft_uncertainty
  expect_equal(u$entropy, log(2), tolerance = 1e-6)
  expect_equal(mt$soft_alignment$f1, 2 / 3, tolerance = 1e-6)

  # ...which is strictly more uncertain than an exact recovery at the
  # same temperature.
  mt2 <- ModuleTester$new(list(r1 = c("a", "b")))
  mt2$set_found(list(f1 = c("a", "b")))
  mt2$align(method = "soft", epsilon = 0.5)
  expect_lt(mt2$soft_uncertainty$entropy, u$entropy)
})

test_that("soft assignment is a valid probabilistic matching", {
  ref <- list(r1 = c("a", "b", "c"), r2 = c("d", "e"), r3 = c("f", "g", "h", "i"))
  found <- list(f1 = c("a", "b", "c"), f2 = c("d", "e", "x"),
                f3 = c("f", "g"), f4 = c("z"))
  mt <- ModuleTester$new(ref)
  mt$set_found(found)
  mt$align(method = "soft", epsilon = 0.1)
  al <- mt$soft_alignment
  P <- mt$soft_assignment
  expect_equal(dim(P), c(3, 4))
  expect_true(all(P >= 0))
  # rows are distributions over found modules (up to dummy mass)
  expect_true(all(rowSums(P) <= 1 + 1e-8))
  # the clearly-best matches carry most mass
  expect_gt(P[1, 1], 0.5)
  expect_gt(P[2, 2], 0.5)
  # expected F1 is between 0 and 1 and below the hard best-match ceiling
  expect_true(all(al$f1 >= 0 & al$f1 <= 1))
  hard <- mt$align(method = "best")
  expect_true(all(al$f1 <= hard$f1 + 1e-9))
  # convergence diagnostics are recorded on the cached entry
  ss <- mt$summary(method = "soft", epsilon = 0.1)
  expect_true(ss$converged)
  expect_true(ss$sinkhorn_iters >= 1L)
})

test_that("soft converges toward hungarian as epsilon decreases", {
  ref <- list(r1 = c("a", "b"), r2 = c("c", "d"), r3 = c("e", "f"))
  found <- list(f1 = c("a", "b", "x"), f2 = c("c", "d"),
                f3 = c("e", "f", "y", "z"))
  mt <- ModuleTester$new(ref)
  mt$set_found(found)
  h <- mt$align(method = "hungarian")
  s <- mt$align(method = "soft", epsilon = 0.05)
  expect_equal(s$f1, h$f1, tolerance = 0.05)
  expect_equal(s$top_found, h$top_found)
})

test_that("invalid solver controls are rejected", {
  mt <- ModuleTester$new(list(r1 = c("a", "b")))
  mt$set_found(list(f1 = c("a", "b")))
  expect_error(mt$align(method = "soft", epsilon = 0), "epsilon")
  expect_error(mt$align(method = "soft", epsilon = -1), "epsilon")
  expect_error(mt$align(method = "soft", epsilon = Inf), "epsilon")
  expect_error(mt$align(method = "soft", epsilon = NA_real_), "epsilon")
  expect_error(mt$summary(method = "soft", epsilon = 0), "epsilon")
  set.seed(1)
  F <- matrix(runif(16), 4, 4)
  expect_error(clr:::.sinkhorn(F, epsilon = 0.5, tol = 0), "tol")
  expect_error(clr:::.sinkhorn(F, epsilon = 0.5, tol = NA_real_), "tol")
  expect_error(clr:::.sinkhorn(F, epsilon = 0.5, max_iter = 0), "max_iter")
  expect_error(clr:::.sinkhorn(F, epsilon = 0.5, max_iter = 10.5), "max_iter")
  expect_error(clr:::.sinkhorn(F, epsilon = 0.5, max_iter = Inf), "max_iter")
})

test_that(".sinkhorn is log-domain stable and reports nonconvergence", {
  set.seed(1)
  F <- matrix(runif(16), 4, 4)

  # converged fit: near-exact uniform margins, diagnostics agree
  sk <- clr:::.sinkhorn(F, epsilon = 0.5)
  expect_true(sk$converged)
  # margins meet the solver's own contract (max deviation < tol = 1e-5)
  expect_lte(max(abs(rowSums(sk$P) - 1 / 4)), 1e-5)
  expect_lte(max(abs(colSums(sk$P) - 1 / 4)), 1e-5)
  expect_true(all(sk$P >= 0))
  expect_lte(sk$marginal_error, 1e-5)

  # tiny epsilon: no NaN/Inf from underflow even without convergence
  sk_tiny <- suppressWarnings(clr:::.sinkhorn(F, epsilon = 1e-3))
  expect_true(all(is.finite(sk_tiny$P)))

  # too few iterations: warns and says so explicitly
  expect_warning(sk_nc <- clr:::.sinkhorn(F, epsilon = 0.5, max_iter = 2L),
                 "did not converge")
  expect_false(sk_nc$converged)
  expect_equal(sk_nc$iters, 2L)
  expect_gt(sk_nc$marginal_error, 0)
})

test_that(".hungarian solves the assignment problem", {
  cost <- matrix(c(4, 1, 3,
                   2, 0, 5,
                   3, 2, 2), 3, 3, byrow = TRUE)
  a <- clr:::.hungarian(cost)
  # optimum: row3->col1 (3), row1->col2 (1), row2->col3... check total
  # brute force over all 6 permutations
  perms <- list(c(1, 2, 3), c(1, 3, 2), c(2, 1, 3),
                c(2, 3, 1), c(3, 1, 2), c(3, 2, 1))
  totals <- vapply(perms, function(p) sum(cost[cbind(p, 1:3)]), numeric(1))
  expect_equal(sum(cost[cbind(a, 1:3)]), min(totals))

  # a second check on 4x4 against exhaustive search
  set.seed(3)
  cost4 <- matrix(runif(16), 4, 4)
  a4 <- clr:::.hungarian(cost4)
  all_perms <- function(n) {
    if (n == 1L) return(list(1L))
    prev <- all_perms(n - 1L)
    out <- list()
    for (p in prev) for (k in seq_len(n)) {
      out[[length(out) + 1L]] <- append(p, n, after = k - 1L)
    }
    out
  }
  totals4 <- vapply(all_perms(4L),
                    function(p) sum(cost4[cbind(p, 1:4)]), numeric(1))
  expect_equal(sum(cost4[cbind(a4, 1:4)]), min(totals4))
})

test_that("headline and soft state coexist; print/summary never clobber", {
  ref <- list(r1 = c("a", "b"), r2 = c("c", "d"))
  mt <- ModuleTester$new(ref)
  mt$set_found(list(f1 = c("a", "b", "x"), f2 = c("c", "d")))

  mt$align(method = "soft", epsilon = 0.5)
  soft_half <- mt$soft_assignment
  # headline alignment is unaffected by the soft run
  expect_equal(mt$alignment$top_found, c("f1", "f2"))
  # a different temperature gets its own cache entry, not a stale reuse
  mt$align(method = "soft", epsilon = 0.05)
  expect_false(isTRUE(all.equal(mt$soft_assignment, soft_half)))
  soft_small <- mt$soft_assignment
  # print() shows the hungarian headline and preserves soft state
  out <- capture.output(mt$print())
  expect_match(out, "hungarian")
  expect_equal(mt$soft_assignment, soft_small)
  # default summary() is the hungarian headline even after soft aligns
  expect_equal(mt$summary()$method, "hungarian")
  # summary(method = "soft") at an uncached temperature computes fresh
  s_new <- mt$summary(method = "soft", epsilon = 0.2)
  expect_equal(s_new$epsilon, 0.2)
  expect_equal(mt$soft_assignment, soft_small)  # binding still at 0.05
})

test_that("find_modules extracts connected components", {
  # two blocks: {1,2,3} connected, {4,5} connected, {6} isolated
  A <- matrix(0, 6, 6)
  A[1, 2] <- A[2, 3] <- A[4, 5] <- 1
  rownames(A) <- paste0("g", 1:6)
  mods <- find_modules(A)
  expect_equal(length(mods), 3)
  expect_equal(sort(mods[[1]]), c("g1", "g2", "g3"))
  expect_equal(sort(mods[[2]]), c("g4", "g5"))
  expect_equal(mods[[3]], "g6")
  # min_size drops the singleton
  mods2 <- find_modules(A, min_size = 2L)
  expect_equal(length(mods2), 2)
  expect_error(find_modules(matrix(1, 2, 3)), "square")
})

test_that("ModuleTester validates its inputs", {
  expect_error(ModuleTester$new(list(c("a", "b"))), "named")
  expect_error(ModuleTester$new(list(r = character(0))), "empty")
  mt <- ModuleTester$new(list(r = c("a", "b")))
  expect_error(mt$align(), "set_found")
  expect_error(mt$set_found(list()), "non-empty")
})

test_that("end-to-end: planted modules recovered from synthetic operator", {
  set.seed(42)
  G <- 12L
  ids <- paste0("g", seq_len(G))
  # plant two modules: genes 1:4 and 5:7 densely linked, rest null
  A <- matrix(0, G, G)
  A[1:4, 1:4] <- 1; A[5:7, 5:7] <- 1
  diag(A) <- 0
  rownames(A) <- colnames(A) <- ids
  mods <- find_modules(A)
  mt <- ModuleTester$new(list(modA = ids[1:4], modB = ids[5:7]))
  mt$set_found(mods)
  s <- mt$summary()
  expect_equal(s$mean_f1, 1)
  expect_equal(s$weighted_f1, 1)
})
