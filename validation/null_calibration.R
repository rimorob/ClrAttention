library(clr); suppressMessages(NULL)
mk <- function(G, n) rbind(matrix(rnorm(G/3*n), G/3), matrix(rt(G/3*n, 3), G/3),
                           matrix(rlnorm(G/3*n), G/3))
res <- NULL
for (tr in c("none", "rank")) for (s in 1:20) {
  set.seed(s); X <- mk(60, 400)
  f <- ClrAttention$new(X)$estimate_mi(threads = 2, transform = tr)$calibrate()
  S <- f$clr_scores; ut <- upper.tri(S)
  suppressMessages(f$select_threshold(B = 20, method = "hc", threads = 2))
  hc <- sum(S[ut] >= f$threshold)
  f$select_threshold(B = 20, method = "fdr", threads = 2)
  res <- rbind(res, data.frame(tr, s, hc, fdr = sum(S[ut] >= f$threshold)))
}
print(aggregate(cbind(hc_any = hc > 0, hc_mean = hc, fdr_any = fdr > 0) ~ tr, res, mean))
