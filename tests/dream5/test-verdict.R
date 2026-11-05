setwd(Sys.getenv("D5_ROOT", unset = "."))
V <- function(d, se) as.character(d5_verdict(d, se)$verdict)
delta <- log2(1.10); z <- qnorm(0.975)

test_that("delta is 10% relative AUPR", expect_equal(D5$delta, log2(1.1)))

test_that("every row of the decision table, with its boundaries", {
  se <- 0.02
  # Superior: p_sup < 0.025 and d >= delta
  expect_equal(V(0.20, se), "Superior")
  expect_equal(V(delta, se), "Superior")                       # d == delta exactly: >= holds
  expect_equal(V(delta - 1e-9, se), "Inconclusive")            # just below delta: significant but below the margin, not equivalent
  # Inferior: p_harm < 0.025 and d <= -delta
  expect_equal(V(-0.20, se), "Inferior")
  expect_equal(V(-delta, se), "Inferior")
  expect_equal(V(-delta + 1e-9, se), "Inconclusive")
  # Equivalent: both non-inferiority tests rejected (95% interval inside +-delta)
  expect_equal(V(0, se), "Equivalent")
  expect_equal(V(0.10, 0.01), "Equivalent")                    # significantly > 0 but inside the margin
  expect_equal(V(-0.10, 0.01), "Equivalent")
  se_edge <- delta / z                                         # interval half-width equal to delta at d = 0
  expect_equal(V(0, se_edge * (1 - 1e-6)), "Equivalent"); expect_equal(V(0, se_edge * (1 + 1e-6)), "Inconclusive")
  expect_equal(V(0.12, 0.02), "Inconclusive")                  # interval crosses +delta
  expect_equal(V(-0.12, 0.02), "Inconclusive")
  # significance boundary of the superiority and harm tests (d well above the margin)
  se_sup <- 0.30 / z
  expect_equal(V(0.30, se_sup * (1 - 1e-6)), "Superior"); expect_equal(V(0.30, se_sup * (1 + 1e-6)), "Inconclusive")
  expect_equal(V(-0.30, se_sup * (1 - 1e-6)), "Inferior"); expect_equal(V(-0.30, se_sup * (1 + 1e-6)), "Inconclusive")
  # large effect but wide interval; zero effect but wide interval
  expect_equal(V(0.30, 0.2), "Inconclusive"); expect_equal(V(0, 0.5), "Inconclusive")
  # d above the margin but the interval reaches below it still counts as Superior (margin applies to d, as written)
  expect_equal(V(0.14, 0.005), "Superior")
  # significant (p_sup < .025) with d between 0 and delta and an interval not inside +-delta: Inconclusive
  expect_equal(V(0.13, 0.01), "Inconclusive")
})

test_that("zero SE and non-finite inputs are handled", {
  expect_equal(V(0.2, 0), "Superior"); expect_equal(V(-0.2, 0), "Inferior"); expect_equal(V(0.05, 0), "Equivalent")
  expect_equal(V(delta, 0), "Superior"); expect_equal(V(-delta, 0), "Inferior"); expect_equal(V(0.2, 0 + 0), "Superior")
  expect_equal(V(NaN, 0.1), "Inconclusive"); expect_equal(V(0.1, NA_real_), "Inconclusive"); expect_equal(V(Inf, 1), "Inconclusive")
})

test_that("p-values follow the formulas of the design", {
  v <- d5_verdict(0.07, 0.05)
  expect_equal(v$p_sup, 1 - pnorm(0.07 / 0.05)); expect_equal(v$p_harm, pnorm(0.07 / 0.05))
  expect_equal(v$p_noninf_lo, 1 - pnorm((0.07 + delta) / 0.05)); expect_equal(v$p_noninf_hi, 1 - pnorm((delta - 0.07) / 0.05))
})

test_that("random (d, SE) cases agree with an independent statement of the table (interval form)", {
  set.seed(42); n <- 20000
  d <- runif(n, -0.4, 0.4); se <- runif(n, 0.005, 0.25)
  ref <- ifelse(d / se > z & d >= delta, "Superior", ifelse(d / se < -z & d <= -delta, "Inferior",
         ifelse((d + delta) / se > z & (delta - d) / se > z, "Equivalent", "Inconclusive")))
  got <- vapply(seq_len(n), function(i) V(d[i], se[i]), "")
  expect_equal(got, ref); expect_setequal(unique(got), c("Superior", "Inferior", "Equivalent", "Inconclusive"))
})
