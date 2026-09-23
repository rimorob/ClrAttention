library(clr)
aupr <- function(s, l) { o <- order(-s); l <- l[o]; tp <- cumsum(l); prec <- tp/seq_along(l); sum(prec[l==1])/sum(l) }
one <- function(seed, tr, bins = "fd") {
  set.seed(seed); n <- 400
  sig <- NULL; mods <- NULL
  for (m in 1:5) { f <- rnorm(n); for (g in 1:6) {
    x <- 0.7*f + sqrt(1-.49)*rnorm(n)
    x <- switch(g %% 3 + 1, x, exp(1.2*x), sign(x)*abs(x)^3)   # monotone distortions
    sig <- rbind(sig, x); mods <- c(mods, m) } }
  nul <- rbind(matrix(rnorm(20*n),20), matrix(rt(20*n,3),20), matrix(rlnorm(20*n),20))
  X <- rbind(sig, nul); G <- nrow(X); cls <- c(rep("sig",30), rep(c("gauss","t3","lnorm"), each=20))
  truth <- outer(c(mods, rep(0,60)), c(mods, rep(-1,60)), "==") * 1
  f <- ClrAttention$new(X)$estimate_mi(threads=2, transform=tr, bins=bins)$calibrate()
  S <- f$clr_scores; ut <- upper.tri(S)
  f$select_threshold(B=20, method="fdr", threads=2); keep <- S >= f$threshold & ut
  fp_cls <- sapply(c("gauss","t3","lnorm"), function(k) sum(keep & !truth & (outer(cls==k, rep(TRUE,G)) | outer(rep(TRUE,G), cls==k))))
  suppressMessages(f$select_threshold(B=20, method="hc", threads=2)); keep_hc <- S >= f$threshold & ut
  c(aupr = aupr(S[ut], truth[ut]), tp_fdr = sum(keep & truth), fp_fdr = sum(keep & !truth),
    fp_fdr_by = fp_cls, tp_hc = sum(keep_hc & truth), fp_hc = sum(keep_hc & !truth), n_true = sum(truth[ut]))
}
for (cfg in list(list("none",10), list("rank",12), list("rank",16), list("rank",20))) {
  r <- sapply(1:6, one, tr = cfg[[1]], bins = cfg[[2]]); cat(cfg[[1]], cfg[[2]], "\n"); print(round(rowMeans(r), 2)) }
