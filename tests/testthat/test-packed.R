set.seed(11)
G <- 90; E <- 70
F <- matrix(rnorm(5 * E), 5)
X <- F[rep(1:5, length.out = G), ] * rep(c(1, 0.6, 0.3), length.out = G) + matrix(rnorm(G * E), G)
rownames(X) <- paste0("g", seq_len(G))
dense <- clr_calibrate(bspline_mi(X, threads = 2), combine = "stouffer")

test_that("packed CLR in double precision is bit-identical to the dense path", {
  bins <- bins_for_genes(X, bins = "hg")
  d2 <- clr:::cpp_clr_packed_double_dense(X, bins, 3L, 1L, 2L)
  expect_identical(unname(d2), unname(dense))
  de <- clr_calibrate(bspline_mi(X, threads = 2), combine = "euclidean")
  expect_identical(unname(clr:::cpp_clr_packed_double_dense(X, bins, 3L, 0L, 2L)), unname(de))
})

test_that("packed float32 CLR matches the dense path closely", {
  p <- clr_packed(X, threads = 2)
  expect_s3_class(p, "clr_packed")
  expect_equal(as.matrix(p), dense, tolerance = 1e-5)
  idx <- c(5, 1, 77, 30)
  expect_equal(clr_packed_block(p, idx), dense[idx, idx], tolerance = 1e-5)
  expect_equal(clr_packed_block(p, rownames(X)[idx]), dense[idx, idx], tolerance = 1e-5)
})

test_that("top-k and the softmax operator reproduce the dense operator", {
  p <- clr_packed(X, threads = 2)
  tk <- clr_packed_topk(p, 20L)
  for (i in c(1, 17, 90)) {
    row <- dense[i, ]; row[i] <- -Inf
    want <- order(row, decreasing = TRUE)[1:20]; want <- want[row[want] > 0]
    expect_identical(tk$idx[i, seq_along(want)], want)
  }
  A_packed <- as.matrix(softmax_operator(p, keff = 10, cap = 50L))
  A_dense <- clr:::.softmax_rows(dense, keff = 10, cap = 50L)
  expect_equal(unname(A_packed), unname(A_dense), tolerance = 1e-5)
  expect_equal(unname(rowSums(A_packed)[rowSums(A_dense) > 0]), rep(1, sum(rowSums(A_dense) > 0)))
})
