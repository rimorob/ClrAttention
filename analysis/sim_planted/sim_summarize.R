#!/usr/bin/env Rscript
# Aggregates results/sim_planted/units/*.rds into docs/plans/sim_planted_results.md following the scoring rules of
# docs/plans/sim_planted_choices.md section 5 (fixed before any run). Usage (repo root): Rscript analysis/sim_planted/sim_summarize.R
args <- commandArgs(TRUE); out_md <- if (length(args)) args[1] else "docs/plans/sim_planted_results.md"
SIZES <- c(5L, 10L, 20L, 40L, 70L, 100L); NS <- if (nzchar(Sys.getenv("SIM_NS"))) as.integer(strsplit(Sys.getenv("SIM_NS"), ",")[[1]]) else c(150L, 400L, 800L)
files <- list.files("results/sim_planted/units", pattern = "rds$", full.names = TRUE)
Us <- lapply(files, readRDS); nunits <- length(Us)
tag <- function(d, u) { d$N <- u$meta$N; d$noise <- u$meta$noise; d$marg <- u$meta$marginal; d$seed <- u$meta$seed; d }
M <- do.call(rbind, lapply(Us, function(u) tag(u$metrics, u)))
D <- do.call(rbind, lapply(Us, function(u) tag(u$depth_metrics, u)))
DM <- do.call(rbind, lapply(Us, function(u) tag(u$depth_mech, u)))
M$id <- paste(M$metric, M$key, M$method, M$N, M$noise, M$marg, sep = "|"); M <- M[order(M$seed), ]
Midx <- split(M$value, M$id)
D$id <- paste(D$metric, D$method, D$t, D$N, D$noise, D$marg, sep = "|"); D <- D[order(D$seed), ]; Didx <- split(D$value, D$id)
vec <- function(metric, key, method, N, noise = "gaussian", marg = "identity") { v <- Midx[[paste(metric, key, method, N, noise, marg, sep = "|")]]; if (is.null(v)) numeric(0) else v }
dvec <- function(metric, method, t, N, noise = "gaussian", marg = "identity") { v <- Didx[[paste(metric, method, t, N, noise, marg, sep = "|")]]; if (is.null(v)) numeric(0) else v }
mechv <- function(arm, field, N, noise = "gaussian", marg = "identity") { us <- Us[order(vapply(Us, function(u) u$meta$seed, 0))]
  us <- us[vapply(us, function(u) u$meta$N == N && u$meta$noise == noise && u$meta$marginal == marg, NA)]
  vapply(us, function(u) u$mech[[arm]][[field]], 0) }
mechmat <- function(arm, field, N, noise = "gaussian", marg = "identity") { us <- Us[order(vapply(Us, function(u) u$meta$seed, 0))]
  us <- us[vapply(us, function(u) u$meta$N == N && u$meta$noise == noise && u$meta$marginal == marg, NA)]
  do.call(rbind, lapply(us, function(u) u$mech[[arm]][[field]])) }

## ---- helpers: scoring rules ---------------------------------------------------------------------------------------
need <- function(n) ceiling(0.9 * n)
gt <- function(a, b) { d <- a - b; length(d) > 0 && mean(d) > 0 && sum(d > 0) >= need(length(d)) }
ge <- function(a, b) !gt(b, a)
apx <- function(a, b) { d <- a - b; length(d) > 0 && abs(mean(d)) <= max(stats::sd(d), 0.10 * max(mean(a), mean(b))) }
between <- function(x, a, b) mean(x) > min(mean(a), mean(b)) && mean(x) < max(mean(a), mean(b))
f3 <- function(x) formatC(as.numeric(x), digits = 3, format = "g")
ms <- function(v) if (!length(v)) "n/a" else sprintf("%s (%s) [%s, %s]", f3(mean(v)), f3(if (length(v) > 1) stats::sd(v) else NA), f3(min(v)), f3(max(v)))
mdtab <- function(df) { df <- as.data.frame(df); df[] <- lapply(df, function(z) gsub("|", "\\|", as.character(z), fixed = TRUE)); c(paste0("| ", paste(names(df), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
  apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))) }

## ---- slopes -------------------------------------------------------------------------------------------------------
slope_rows <- list()
for (u in Us) { m <- u$metrics[u$metrics$metric == "mod_mean_score", ]; sz <- u$meta$sizes; lam <- u$meta$lam
  for (me in unique(m$method)) { mm <- m[m$method == me, ]; y <- mm$value[order(as.integer(mm$key))]
    for (fit in c("n20", "all")) { ix <- if (fit == "n20") sz >= 20 else rep(TRUE, length(sz))
      for (sc in c("log", "raw")) { yy <- if (sc == "log") log(y[ix]) else y[ix]
        sl <- if (all(is.finite(yy))) unname(stats::coef(stats::lm(yy ~ log(sz[ix]) + factor(lam[ix])))[2]) else NA_real_
        slope_rows[[length(slope_rows) + 1L]] <- data.frame(N = u$meta$N, noise = u$meta$noise, marg = u$meta$marginal, seed = u$meta$seed, method = me, fit = fit, scale = sc, slope = sl) } } } }
SL <- do.call(rbind, slope_rows); SL <- SL[order(SL$seed), ]; SLidx <- split(SL$slope, paste(SL$method, SL$fit, SL$scale, SL$N, SL$noise, SL$marg, sep = "|"))
slv <- function(method, fit, scale, N, noise = "gaussian", marg = "identity") { v <- SLidx[[paste(method, fit, scale, N, noise, marg, sep = "|")]]; if (is.null(v)) numeric(0) else v }

## ---- prediction bookkeeping ---------------------------------------------------------------------------------------
CL <- list(); QT <- list()
clause <- function(pred, label, N, ok) CL[[length(CL) + 1L]] <<- data.frame(pred = pred, clause = label, N = N, ok = isTRUE(ok), stringsAsFactors = FALSE)
quant <- function(pred, label, N, v) QT[[length(QT) + 1L]] <<- data.frame(pred = pred, quantity = label, N = N, cell = ms(v), stringsAsFactors = FALSE)
verdict_of <- function(ok) if (!length(ok)) "not evaluable" else if (all(ok)) "held" else if (!any(ok)) "did not hold" else "partly held"
T2 <- function(me, b, N, ...) vec("aupr_T2", b, me, N, ...)
BIG <- c(20L, 40L, 70L, 100L)
growth_ok <- function(g) { inv <- sum(outer(seq_along(g), seq_along(g), function(i, j) i < j & g[j] < g[i])); g[length(g)] > g[1] && inv <= 1 }

for (N in NS) {
  ## R1
  gap_cm <- gap_tm <- numeric(0)
  for (b in BIG) { tc <- T2("tcaclr", b, N); cl <- T2("clr", b, N); ac <- T2("aclr", b, N)
    clause("R1", sprintf("T2 n=%d: tcaCLR > CLR", b), N, gt(tc, cl)); clause("R1", sprintf("T2 n=%d: CLR > aCLR", b), N, gt(cl, ac))
    gap_cm <- c(gap_cm, mean(cl - ac)); gap_tm <- c(gap_tm, mean(tc - cl))
    for (me in c("tcaclr", "clr", "aclr")) quant("R1", sprintf("T2 AUPR %s n=%d", me, b), N, T2(me, b, N)) }
  clause("R1", "aCLR deficit (CLR - aCLR) grows with n", N, growth_ok(gap_cm)); clause("R1", "tcaCLR - CLR gap grows with n", N, growth_ok(gap_tm))
  sa <- slv("aclr", "n20", "log", N); sc <- slv("clr", "n20", "log", N); st <- slv("tcaclr", "n20", "raw", N)
  clause("R1", "slope aCLR in [-1.25, -0.75]", N, mean(sa) >= -1.25 && mean(sa) <= -0.75); clause("R1", "slope CLR in [-0.75, -0.25]", N, mean(sc) >= -0.75 && mean(sc) <= -0.25)
  clause("R1", "|slope| tcaCLR < 0.3 z per nat", N, abs(mean(st)) < 0.3)
  quant("R1", "slope aCLR (log score on log n, n>=20)", N, sa); quant("R1", "slope CLR (log score on log n, n>=20)", N, sc); quant("R1", "slope tcaCLR (z per nat, n>=20)", N, st)
  sc_raw <- slv("clr", "n20", "raw", N); st_log <- slv("tcaclr", "n20", "log", N)
  quant("R1", "memo failure test, z units: |slope tcaCLR| - |slope CLR| (negative = tcaCLR flatter)", N, abs(st) - abs(sc_raw))
  quant("R1", "memo failure test, log-score units: |slope tcaCLR| - |slope CLR| (negative = tcaCLR flatter)", N, abs(st_log) - abs(sc))
  ## R2
  for (me in c("tcaclr", "clr", "aclr")) { clause("R2", sprintf("%s: bin 5 below bin 20", me), N, gt(T2(me, 20L, N), T2(me, 5L, N)))
    clause("R2", sprintf("%s: bin 10 intermediate", me), N, between(T2(me, 10L, N), T2(me, 5L, N), T2(me, 20L, N)))
    for (b in c(5L, 10L)) quant("R2", sprintf("T2 AUPR %s n=%d", me, b), N, T2(me, b, N)) }
  clause("R2", "n=5: CLR >= tcaCLR", N, ge(T2("clr", 5L, N), T2("tcaclr", 5L, N))); clause("R2", "n=5: tcaCLR ~ aCLR", N, apx(T2("tcaclr", 5L, N), T2("aclr", 5L, N)))
  ## R3
  clause("R3", "T3 pooled: aCLR ~ tcaCLR", N, apx(vec("aupr_T3", "all", "aclr", N), vec("aupr_T3", "all", "tcaclr", N)))
  clause("R3", "T3 pooled: tcaCLR > CLR", N, gt(vec("aupr_T3", "all", "tcaclr", N), vec("aupr_T3", "all", "clr", N)))
  clause("R3", "T3 pooled: aCLR > CLR", N, gt(vec("aupr_T3", "all", "aclr", N), vec("aupr_T3", "all", "clr", N)))
  for (me in c("aclr", "tcaclr", "clr")) quant("R3", sprintf("T3 AUPR pooled %s", me), N, vec("aupr_T3", "all", me, N))
  skeys <- unique(M$key[M$metric == "aupr_T3s"]); gaps <- vapply(skeys, function(k) mean((vec("aupr_T3s", k, "aclr", N) + vec("aupr_T3s", k, "tcaclr", N)) / 2 - vec("aupr_T3s", k, "clr", N)), 0)
  clause("R3", "largest gap (mean(aCLR, tcaCLR) - CLR) at m=40, b=0.35", N, names(gaps)[which.max(gaps)] == "m40_b0.35")
  for (k in skeys) quant("R3", sprintf("stratum %s: mean(aCLR,tcaCLR) - CLR (seed-mean AUPR gap)", k), N, (vec("aupr_T3s", k, "aclr", N) + vec("aupr_T3s", k, "tcaclr", N)) / 2 - vec("aupr_T3s", k, "clr", N))
  ## R4
  a4 <- vec("aupr_T4", "all", "aclr", N); t4 <- vec("aupr_T4", "all", "tcaclr", N); c4 <- vec("aupr_T4", "all", "clr", N); r4 <- vec("aupr_T4", "all", "absr", N)
  clause("R4", "T4: CLR >= tcaCLR", N, ge(c4, t4)); clause("R4", "T4: tcaCLR >= aCLR", N, ge(t4, a4))
  clause("R4", "T4: |r| >= 0.8 x best of CLR, aCLR, tcaCLR", N, mean(r4) >= 0.8 * max(mean(c4), mean(a4), mean(t4)))
  for (nm in c("clr", "tcaclr", "aclr", "absr")) quant("R4", sprintf("T4 AUPR %s", nm), N, vec("aupr_T4", "all", nm, N))
  ## R5
  slc_g <- slv("clr", "n20", "log", N); slc_t <- slv("clr", "n20", "log", N, "t4")
  gapf <- function(noise) mean(vapply(BIG, function(b) mean(T2("tcaclr", b, N, noise) - T2("clr", b, N, noise)), 0))
  clause("R5", "heavy tail: seed SD of CLR slope larger than Gaussian", N, length(slc_t) > 1 && stats::sd(slc_t) > stats::sd(slc_g))
  clause("R5", "heavy tail: tcaCLR - CLR gap (T2, n>=20, mean over bins) larger than Gaussian", N, gapf("t4") > gapf("gaussian"))
  quant("R5", "CLR slope, Gaussian noise", N, slc_g); quant("R5", "CLR slope, Student-t noise", N, slc_t)
  quant("R5", "tcaCLR - CLR T2 gap n>=20 (per seed, mean over bins), Gaussian", N, rowMeans(sapply(BIG, function(b) T2("tcaclr", b, N) - T2("clr", b, N))))
  quant("R5", "tcaCLR - CLR T2 gap n>=20 (per seed, mean over bins), Student-t", N, rowMeans(sapply(BIG, function(b) T2("tcaclr", b, N, "t4") - T2("clr", b, N, "t4"))))
  dI <- mean(vec("aupr_T1", "all", "tcaclr", N) - vec("aupr_T1", "all", "clr", N)); dC <- mean(vec("aupr_T1", "all", "tcaclr", N, marg = "copula") - vec("aupr_T1", "all", "clr", N, marg = "copula"))
  clause("R5", "copula: sign of tcaCLR - CLR on T1 equals identity arm", N, sign(dI) == sign(dC))
  clause("R5", "copula: |r| T1 below CLR", N, gt(vec("aupr_T1", "all", "clr", N, marg = "copula"), vec("aupr_T1", "all", "absr", N, marg = "copula")))
  clause("R5", "copula: |r| T1 below tcaCLR", N, gt(vec("aupr_T1", "all", "tcaclr", N, marg = "copula"), vec("aupr_T1", "all", "absr", N, marg = "copula")))
  for (nm in c("absr", "clr", "aclr", "tcaclr")) { quant("R5", sprintf("T1 AUPR %s, identity", nm), N, vec("aupr_T1", "all", nm, N)); quant("R5", sprintf("T1 AUPR %s, copula", nm), N, vec("aupr_T1", "all", nm, N, marg = "copula")) }
  ## R7
  for (b in BIG) { clause("R7", sprintf("T2 n=%d: sinkhorn raw ~ aCLR", b), N, apx(T2("sink_raw", b, N), T2("aclr", b, N))); quant("R7", sprintf("T2 AUPR sinkhorn raw n=%d", b), N, T2("sink_raw", b, N)); quant("R7", sprintf("T2 AUPR aCLR n=%d", b), N, T2("aclr", b, N)) }
  clause("R7", "hub share of top 1%: sinkhorn raw < aCLR", N, gt(vec("hubshare_pi", "top1", "aclr", N), vec("hubshare_pi", "top1", "sink_raw", N)))
  quant("R7", "hub share top 1% (high-pi set), aCLR", N, vec("hubshare_pi", "top1", "aclr", N)); quant("R7", "hub share top 1% (high-pi set), sinkhorn raw", N, vec("hubshare_pi", "top1", "sink_raw", N))
  ## R8
  hsT1 <- vec("aupr_T1", "all", "hubscore", N); quant("R8", "T1 AUPR of the hub-score reference (pi_i + pi_j)", N, hsT1)
  rises <- list()
  for (me in c("aclr", "tcaclr")) { a9 <- dvec("aupr_T1", me, 9, N); a36 <- dvec("aupr_T1", me, 36, N); a72 <- dvec("aupr_T1", me, 72, N)
    h9 <- dvec("hubshare_pi", me, 9, N); h36 <- dvec("hubshare_pi", me, 36, N); h72 <- dvec("hubshare_pi", me, 72, N)
    clause("R8", sprintf("%s: T1 at t=36 below t=9", me), N, gt(a9, a36)); clause("R8", sprintf("%s: T1 at t=72 below t=9", me), N, gt(a9, a72))
    clause("R8", sprintf("%s: T1(72) closer to hub-score level than T1(9)", me), N, gt(abs(log(a9 / hsT1)), abs(log(a72 / hsT1))))
    clause("R8", sprintf("%s: hub share at t=36 above t=9", me), N, gt(h36, h9)); clause("R8", sprintf("%s: hub share at t=72 above t=9", me), N, gt(h72, h9))
    rises[[me]] <- h72 - h9
    for (t in c(9, 36, 72)) { quant("R8", sprintf("%s T1 AUPR t=%d", me, t), N, dvec("aupr_T1", me, t, N)); quant("R8", sprintf("%s hub share top 1%% t=%d", me, t), N, dvec("hubshare_pi", me, t, N)) } }
  clause("R8", "hub-share rise 9 -> 72: aCLR larger than tcaCLR", N, gt(rises$aclr, rises$tcaclr))
  ## R9
  for (b in BIG) { ks <- which(Us[[1]]$meta$sizes == b)
    sa_v <- function(me) { x <- M[M$metric == "sat_auc" & M$method == me & M$N == N & M$noise == "gaussian" & M$marg == "identity" & suppressWarnings(as.integer(M$key)) %in% ks, ]; tapply(x$value, x$seed, mean) }
    va <- sa_v("aclr"); vc <- sa_v("clr"); vt <- sa_v("tcaclr")
    clause("R9", sprintf("n=%d: aCLR ranks satellite-member above member-member (P > 0.5)", b), N, mean(va) > 0.5 && sum(va > 0.5) >= need(length(va)))
    clause("R9", sprintf("n=%d: CLR ranks member-member higher (P < 0.5)", b), N, mean(vc) < 0.5 && sum(vc < 0.5) >= need(length(vc)))
    clause("R9", sprintf("n=%d: tcaCLR ranks member-member higher (P < 0.5)", b), N, mean(vt) < 0.5 && sum(vt < 0.5) >= need(length(vt)))
    quant("R9", sprintf("P(satellite-member > member-member), n=%d, aCLR", b), N, va); quant("R9", sprintf("P(...), n=%d, CLR", b), N, vc); quant("R9", sprintf("P(...), n=%d, tcaCLR", b), N, vt) }
  ## R10
  s10 <- mechv("t9", "sigma_median", N); s40 <- mechv("k40", "sigma_median", N)
  g10 <- vec("aupr_T1", "all", "tcaclr", N) - vec("aupr_T1", "all", "aclr", N); g40 <- vec("aupr_T1", "all", "tcaclr_k40", N) - vec("aupr_T1", "all", "aclr_k40", N)
  clause("R10", "median sigma_bg at keff 40 below keff 10", N, gt(s10, s40)); clause("R10", "tcaCLR - aCLR T1 gap smaller at keff 40 than at keff 10", N, gt(g10, g40))
  quant("R10", "median sigma_bg (nats), keff 10", N, s10); quant("R10", "median sigma_bg (nats), keff 40", N, s40); quant("R10", "T1 gap tcaCLR - aCLR, keff 10", N, g10); quant("R10", "T1 gap tcaCLR - aCLR, keff 40", N, g40)
  ## R11
  rel <- function(me) { v <- vec("aupr_T1", "all", me, N); stats::sd(v) / mean(v) }
  for (me in c("clr", "aclr", "tcaclr")) clause("R11", sprintf("relative SD of T1 AUPR: pi log-lift larger than %s", me), N, rel("pi_lift") > rel(me))
  for (me in c("aclr", "tcaclr")) clause("R11", sprintf("hub share of top 1%%: pi log-lift above %s", me), N, gt(vec("hubshare_pi", "top1", "pi_lift", N), vec("hubshare_pi", "top1", me, N)))
  clause("R11", "column-mean log-lift T1 strictly between aCLR and tcaCLR", N, between(vec("aupr_T1", "all", "colmean_lift", N), vec("aupr_T1", "all", "aclr", N), vec("aupr_T1", "all", "tcaclr", N)))
  for (me in c("pi_lift", "colmean_lift", "aclr", "tcaclr", "clr")) { quant("R11", sprintf("T1 AUPR %s", me), N, vec("aupr_T1", "all", me, N)); quant("R11", sprintf("relative SD of T1 AUPR %s (single value across seeds)", me), N, rel(me)) }
  quant("R11", "hub share top 1%, pi log-lift", N, vec("hubshare_pi", "top1", "pi_lift", N))
  ## R12
  for (me in c("infl_1.5", "infl_2")) { clause("R12", sprintf("T1: tcaCLR > %s (inflation does not reach tcaCLR)", me), N, gt(vec("aupr_T1", "all", "tcaclr", N), vec("aupr_T1", "all", me, N))); quant("R12", sprintf("T1 AUPR %s", me), N, vec("aupr_T1", "all", me, N)) }
  clause("R12", "hub-stratified T1 (true hub loading sum): sinkhorn -> tcaCLR >= tcaCLR", N, ge(vec("aupr_T1_hub_true", "all", "sink_tca", N), vec("aupr_T1_hub_true", "all", "tcaclr", N)))
  for (me in c("sink_tca", "tcaclr")) { quant("R12", sprintf("hub-stratified T1 (true hub loading sum) %s", me), N, vec("aupr_T1_hub_true", "all", me, N)); quant("R12", sprintf("hub-stratified T1 (MI row mean) %s [reported]", me), N, vec("aupr_T1_hub_mi", "all", me, N)) }
  quant("R12", "T1 AUPR tcaCLR", N, vec("aupr_T1", "all", "tcaclr", N))
}
## R6
mx <- max(vapply(Us, function(u) max(u$ident$maxdiff_Z, u$ident$maxdiff_S), 0)); rle <- vapply(Us, function(u) u$ident$rows_le10, 0)
clause("R6", "max |difference| <= 1e-10 in every compendium", 0L, mx <= 1e-10)

CLd <- do.call(rbind, CL); QTd <- do.call(rbind, QT)
Rs <- paste0("R", 1:12)
pred_text <- c(
 R1 = "(T2, n >= 20): tcaCLR > CLR > aCLR, with the aCLR deficit and the tcaCLR-CLR gap both growing with n; slopes ~ -1 (aCLR), ~ -0.5 (CLR), |slope| < 0.3 z/nat (tcaCLR).",
 R2 = "(T2, n = 5): CLR >= tcaCLR ~ aCLR; all three below their n = 20 level (forced leakage); n = 10 intermediate.",
 R3 = "(T3, loose regulons): aCLR ~ tcaCLR > CLR, gap largest for b = 0.35 and m = 40.",
 R4 = "(T4, direct edges): CLR >= tcaCLR >= aCLR; |r| competitive.",
 R5 = "(heavy-tail arm): CLR's size penalty becomes noisier (sigma_i inflated), tcaCLR's advantage on n >= 20 grows; (copula arm): CLR and tcaCLR unchanged in rank, |r| loses.",
 R6 = "(identity): tcaCLR at t = 1 equals Stouffer of within-candidate row z of CLR to 1e-10.",
 R7 = "(lemma): Sinkhorn-balanced walk, raw readout, ties aCLR on n >= 20 (within seed SD) and reduces the hub share of the top 1%.",
 R8 = "(depth): at t = 36 and 72 all methods' T1 falls toward the hub-score level and the hub share of the top 1% rises; aCLR's rises faster than tcaCLR's.",
 R9 = "(satellites): aCLR ranks satellite-member pairs above member-member pairs of the same module for n >= 20; CLR and tcaCLR rank member-member higher.",
 R10 = "(wide-row arm, keff 40): sigma_bg falls, tcaCLR's T1 gain over aCLR shrinks or vanishes; CLR unaffected.",
 R11 = "(lift): pi log-lift unstable (large seed SD, hub-dominated top 1%); column-mean log-lift between aCLR and tcaCLR on T1.",
 R12 = "(inflation control): fixed-exponent inflation does not reach tcaCLR on T1 for both r; Sinkhorn(A) -> walk -> tcaCLR is at or above tcaCLR on hub-stratified T1 (speculation; a failure here is informative).")

## ---- H2 bin-edge rule ---------------------------------------------------------------------------------------------
h2 <- do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(SIZES, function(b) { d <- T2("tcaclr", b, N) - T2("aclr", b, N)
  data.frame(N = N, n = b, positive_seeds = sum(d > 0), seeds = length(d), mean_diff = mean(d), qual = sum(d > 0) >= need(length(d))) }))))
edge <- vapply(NS, function(N) { x <- h2[h2$N == N & h2$qual, ]; if (nrow(x)) min(x$n) else NA_real_ }, 0)
contig <- vapply(NS, function(N) { e <- edge[match(N, NS)]; if (is.na(e)) NA else all(h2$qual[h2$N == N & h2$n >= e]) }, NA)

## ---- write --------------------------------------------------------------------------------------------------------
L <- character(0); w <- function(...) L <<- c(L, paste0(...))
w("# Planted-module simulation: results"); w("")
w("Specification: `docs/plans/fable_theory_memo.md` section 5.4. Implementation choices and the scoring rules for R1 to R12 were fixed and committed before any run (`docs/plans/sim_planted_choices.md`); the predictions below are quoted unchanged from the memo. Synthetic data only. No method parameter was tuned (keff 10, cap 50, alpha 0.5, t = 9; the keff 40 / cap 200 arm and the depth grid are descriptive). Code: `analysis/sim_planted/`; per-compendium results (not committed): `results/sim_planted/units/`.")
w(""); w(sprintf("Compendia analysed: %d of 120. Values are mean (SD across seeds) [min, max] over the 10 seeds of the cell; the primary arm is Gaussian noise with identity marginals, shown at N = 150, 400, 800. Compendia per cell:", nunits))
cells <- as.data.frame(table(N = vapply(Us, function(u) u$meta$N, 0), noise = vapply(Us, function(u) u$meta$noise, ""), marginal = vapply(Us, function(u) u$meta$marginal, "")))
cells <- cells[cells$Freq > 0, ]; names(cells)[4] <- "compendia"; w(""); L <- c(L, mdtab(cells)); w("")
w("## Verdicts"); w(""); w("A prediction is judged on clause cells (a clause evaluated at one N); held = every cell true, did not hold = none true, partly held = otherwise (rules in the choices file, section 5).")
w("")
vt <- do.call(rbind, lapply(Rs, function(r) { x <- CLd[CLd$pred == r, ]; data.frame(Prediction = r, Verdict = verdict_of(x$ok), Cells_true = sprintf("%d of %d", sum(x$ok), nrow(x)), Statement = pred_text[[r]]) }))
L <- c(L, mdtab(vt)); w("")
w("## H2 bin-edge rule for `docs/plans/tcaclr_confirm_plan.md`"); w("")
w("Rule: in the primary arm (Gaussian, identity), the smallest module size at which tcaCLR minus aCLR (T2 per-size-bin pooled AUPR) is positive in at least 9 of the 10 seeds, evaluated at each N. Strictly positive difference per seed.")
w(""); L <- c(L, mdtab(data.frame(N = h2$N, `module size` = h2$n, `seeds with tcaCLR > aCLR` = sprintf("%d of %d", h2$positive_seeds, h2$seeds), `mean difference in AUPR` = f3(h2$mean_diff), `meets 9 of 10` = ifelse(h2$qual, "yes", "no"), check.names = FALSE))); w("")
w("Edge per N (smallest qualifying size): ", paste(sprintf("N = %d: %s", NS, ifelse(is.na(edge), "not reached", as.character(edge))), collapse = "; "), ". Every larger size also qualifies: ", paste(sprintf("N = %d: %s", NS, ifelse(is.na(contig), "n/a", ifelse(contig, "yes", "no"))), collapse = "; "), ".")
ov <- if (any(is.na(edge))) "not reached at every N" else as.character(max(edge))
w(""); w("**Overall (the smallest size that qualifies at every N): ", ov, "**; the most favourable N gives ", if (all(is.na(edge))) "none" else as.character(min(edge, na.rm = TRUE)), ".")
w("")
nf <- "analysis/sim_planted/results_notes.md"; if (file.exists(nf)) { w("## Plain statements on what held and what did not"); w(""); L <- c(L, readLines(nf)); w("") }
for (r in Rs) {
  w("## ", r); w(""); w("Prediction (memo): ", pred_text[[r]]); w("")
  x <- CLd[CLd$pred == r, ]; w("Verdict: **", verdict_of(x$ok), "** (", sum(x$ok), " of ", nrow(x), " clause cells true).")
  w(""); q <- QTd[QTd$pred == r, ]
  if (nrow(q)) { qw <- reshape(q[, c("quantity", "N", "cell")], idvar = "quantity", timevar = "N", direction = "wide"); names(qw) <- sub("cell.", "N = ", names(qw)); qw <- qw[match(unique(q$quantity), qw$quantity), ]; L <- c(L, mdtab(qw)); w("") }
  cw <- reshape(x[, c("clause", "N", "ok")], idvar = "clause", timevar = "N", direction = "wide"); names(cw) <- sub("ok.", "N = ", names(cw))
  cw[-1] <- lapply(cw[-1], function(z) ifelse(z, "true", "false")); cw <- cw[match(unique(x$clause), cw$clause), ]; L <- c(L, mdtab(cw)); w("")
  fl <- x[!x$ok, ]
  if (nrow(fl)) { w("Failed clauses: ", paste(sprintf("%s (N = %s)", fl$clause, ifelse(fl$N == 0, "all", fl$N)), collapse = "; "), "."); w("") }
  if (r == "R6") { w(sprintf("Maximum absolute difference over %d compendia: %.2e; rows with at most 10 candidates (excluded): %s per compendium (range).", nunits, mx, paste(range(rle), collapse = " to "))); w("") }
}
## mechanism curves
w("## Mechanism curves (primary arm)"); w("")
w("Mean within-module pair score by module size (seed means; aCLR and CLR on the raw scale, tcaCLR in z units), with fitted slopes (module-level regression on log n with a lambda fixed effect; n >= 20 primary, all sizes secondary).")
for (N in NS) { w(""); w("N = ", N); w("")
  tb <- do.call(rbind, lapply(c("aclr", "clr", "tcaclr", "absr"), function(me) { r <- sapply(SIZES, function(b) { ks <- which(Us[[1]]$meta$sizes == b); x <- M[M$metric == "mod_mean_score" & M$method == me & M$N == N & M$noise == "gaussian" & M$marg == "identity" & suppressWarnings(as.integer(M$key)) %in% ks, ]; mean(x$value) })
    sl <- function(fit, sc) ms(slv(me, fit, sc, N)); data.frame(method = me, t(setNames(f3(r), paste0("n=", SIZES))), `slope n>=20 (log)` = sl("n20", "log"), `slope n>=20 (raw)` = sl("n20", "raw"), `slope all (log)` = sl("all", "log"), check.names = FALSE) }))
  L <- c(L, mdtab(tb)) }
w(""); w("## Per-row sigma_bg, retention m(n), hub share versus t (primary arm)"); w("")
w("Median per-row sigma_bg (SD of log P^t over reachable off-diagonal entries, nats), coefficient of variation across rows, and unreachable share, by depth t (tcaCLR operator, keff 10):"); w("")
sb <- do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(c(3, 6, 9, 18, 36, 72), function(t) { x <- DM[DM$N == N & DM$noise == "gaussian" & DM$marg == "identity" & DM$t == t, ]
  data.frame(N = N, t = t, `median sigma_bg` = ms(x$sigma_median), `CV across rows` = ms(x$cv), `unreachable share` = ms(x$unreach_share), check.names = FALSE) }))))
L <- c(L, mdtab(sb)); w("")
w("sigma_bg at t = 1 and t = 9 for the keff 10, keff 40 and Sinkhorn operators (median across rows; CV across rows):"); w("")
sg <- do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(c("t1", "t9", "k40", "sink"), function(a) data.frame(N = N, operator = c(t1 = "keff 10, t = 1", t9 = "keff 10, t = 9", k40 = "keff 40, t = 9", sink = "Sinkhorn, t = 9")[[a]], `median sigma_bg` = ms(mechv(a, "sigma_median", N)), `CV` = ms(mechv(a, "sigma_cv", N)), `unreachable share` = ms(mechv(a, "unreach_share", N)), check.names = FALSE)))))
L <- c(L, mdtab(sg)); w("")
w("In-module retention by module size: the share of attention A on the other members (t = 1, A only) and the mass sum_{j in module} P^t_ij (self included) at t = 1 and t = 9 (mean over members and modules of that size, then seeds):"); w("")
rt <- do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(SIZES, function(b) { ks <- which(Us[[1]]$meta$sizes == b)
  g <- function(arm, fld) ms(rowMeans(mechmat(arm, fld, N)[, ks, drop = FALSE]))
  data.frame(N = N, n = b, `A in-module share` = g("t1", "A_inshare"), `m_1(n)` = g("t1", "retention"), `m_9(n)` = g("t9", "retention"), `m_9(n), keff 40` = g("k40", "retention"), check.names = FALSE) }))))
L <- c(L, mdtab(rt)); w("")
w("Hub share of the top 1% of pairs (pairs containing at least one of the 30 highest-pi genes) and T1 AUPR versus depth; chance level 2.0%, hub-score reference in R8:"); w("")
hd <- do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(c("aclr", "tcaclr"), function(me) do.call(rbind, lapply(c(3, 6, 9, 18, 36, 72), function(t)
  data.frame(N = N, method = me, t = t, `T1 AUPR` = ms(dvec("aupr_T1", me, t, N)), `hub share top 1%` = ms(dvec("hubshare_pi", me, t, N)), `planted-hub share` = ms(dvec("hubshare_planted", me, t, N)), check.names = FALSE)))))))
L <- c(L, mdtab(hd)); w("")
w("## All methods, primary arm (pooled AUPR; base rates T1 ", f3(mean(vec("base_T1", "all", "clr", 400))), ", T3 ", f3(mean(vec("base_T3", "all", "clr", 400))), ", T4 ", f3(mean(vec("base_T4", "all", "clr", 400))), ")"); w("")
meths <- c("absr", "clr", "aclr", "tcaclr", "tcaclr_t1", "sink_raw", "sink_tca", "pi_lift", "colmean_lift", "infl_1.5", "infl_2", "aclr_k40", "tcaclr_k40")
for (N in NS) { w("N = ", N); w(""); L <- c(L, mdtab(do.call(rbind, lapply(meths, function(me) data.frame(method = me, T1 = ms(vec("aupr_T1", "all", me, N)), T3 = ms(vec("aupr_T3", "all", me, N)), T4 = ms(vec("aupr_T4", "all", me, N)),
  `T1 hub-strat (true)` = ms(vec("aupr_T1_hub_true", "all", me, N)), `T1 hub-strat (MI)` = ms(vec("aupr_T1_hub_mi", "all", me, N)), `satellite share top 1%` = ms(vec("sat_share_top1", "top1", me, N)), check.names = FALSE))))); w("") }
w("T2 per-bin AUPR (primary arm, seed means; per-seed values in the unit files):"); w("")
L <- c(L, mdtab(do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(c("absr", "clr", "aclr", "tcaclr", "sink_raw", "aclr_k40", "tcaclr_k40"), function(me) data.frame(N = N, method = me, t(setNames(f3(sapply(SIZES, function(b) mean(T2(me, b, N)))), paste0("n=", SIZES))), check.names = FALSE))))))); w("")
w("T2 per-bin AUPR in the other arms (tcaCLR minus aCLR and CLR minus aCLR, seed mean):"); w("")
arms <- expand.grid(noise = c("gaussian", "t4"), marg = c("identity", "copula"), stringsAsFactors = FALSE)
L <- c(L, mdtab(do.call(rbind, lapply(NS, function(N) do.call(rbind, lapply(seq_len(nrow(arms)), function(i) { a <- arms[i, ]
  do.call(rbind, lapply(c("tcaclr", "clr"), function(me) data.frame(N = N, arm = paste(a$noise, a$marg), contrast = paste(me, "- aclr"), t(setNames(f3(sapply(SIZES, function(b) mean(T2(me, b, N, a$noise, a$marg) - T2("aclr", b, N, a$noise, a$marg)))), paste0("n=", SIZES))), check.names = FALSE))) }))))))
w(""); w("## Run record"); w("")
w(sprintf("Compendia %d; wall per compendium %s s (R process peak %s GB).", nunits, ms(vapply(Us, function(u) u$meta$wall_s, 0)), ms(vapply(Us, function(u) u$meta$vmhwm_gb, 0))))
sk <- vapply(Us, function(u) u$mech$sink$sinkhorn_resid, 0); w(sprintf("Sinkhorn balancing of the lazy operator hit the 2,000-iteration cap in %d of %d compendia; final max |column sum - 1| %s (tolerance 1e-9 not reached where the cap was hit; the walk used the last iterate).", sum(vapply(Us, function(u) u$mech$sink$sinkhorn_iters, 0) >= 2000), nunits, ms(sk)))
writeLines(L, out_md); cat(sprintf("wrote %s (%d lines), %d units\n", out_md, length(L), nunits))
