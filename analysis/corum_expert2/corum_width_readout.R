#!/usr/bin/env Rscript
# docs/plans/corum_expert2.md. EXPLORATORY PATHFINDING ON AN IN-SAMPLE TRUTH, NOT PREREGISTERED,
# NOTHING CONFIRMATORY. Nothing from this script goes into paper/ or any paper table.
#
# K562 CORUM was used in analysis/corum_exploratory/, the keff test, the CLR(log P^t) test and
# analysis/rclr_k562/, so it is not held out. Same data as those (half A, 5,671 truth genes x
# 707 knockdown profiles), same CLR (B-spline MI, 'hg' bins, Stouffer), same lazy walk
# (alpha 0.5), same symmetric readout S = (P^t + P^t')/2 unless a readout is named.
#
# Task A questions (operator and readout lists fixed before running):
#   (1) cap and width beyond the keff test's grid (keff <= 40, cap 50 fixed): keff 10 at cap 100
#       and 200; keff 25/40 at cap 50; keff 40 at cap 100; keff 80 at cap 100 and 200; uniform
#       attention over the top 50 / 100 positive CLR partners (a kNN walk).
#   (2) readouts of the reference operator (K10, cap 50) that are less scale-dependent or
#       one-sidedness-proof: min(P, P'); mean of row percentiles; degree-scaled mass
#       (d_i P_ij + d_j P_ji)/2 with d_i the number of CLR partners passing a permutation BH
#       q = 0.05 threshold (label-free); a3r (reference, docs/plans/corum_analysis.md).
#       Control: row-percentile readout of CLR itself.
#   (3) evaluation: raw CORUM, overlap-merged CORUM (merge_sets of analysis/clr_logpt_test.R),
#       merged CORUM without RP-dominated complexes, small (3-6) and large (>= 20) raw complexes,
#       per-complex (macro) AP by size bin, hub-stratified pooled AUPR; guards DoRothEA A-B and
#       CollecTRI (also in-sample K562 truths).
# Task B (label-free): row diagnostics of MI and of |r| against a within-gene permutation null.
#
# Usage (repo root): Rscript analysis/corum_expert2/corum_width_readout.R [--out results/corum_expert2/k562] [--threads 8] [--jack 100] [--quick 0]
args <- commandArgs(TRUE)
opt <- list(out = "results/corum_expert2/k562", threads = 8L, jack = 100L, quick = 0L)
for (k in seq(1, length(args), by = 2)) { key <- sub("^--", "", args[k]); opt[[key]] <- if (key == "out") args[k + 1] else as.integer(args[k + 1]) }
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
source("analysis/corum_expert2/null_diag.R")
QUICK <- opt$quick > 0L

## ---- data ----------------------------------------------------------------------------
raw <- utils::read.delim("results/baselines_k562ess_truthgenes_genie3/expression.tsv", check.names = FALSE, stringsAsFactors = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- raw[[1]]; storage.mode(X) <- "double"; rm(raw)
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]
if (QUICK) { v <- apply(X, 1, stats::var); X <- X[sort(order(-v)[seq_len(opt$quick)]), ]; rm(v) }
genes <- rownames(X); G <- length(genes)
say(sprintf("%d genes x %d profiles (quick = %d)", G, ncol(X), opt$quick))
fit_of <- function(Xm) { f <- ClrAttention$new(Xm)$estimate_mi(bins = "hg", transform = "none", threads = opt$threads)
  f$calibrate(method = "normal", combine = "stouffer"); f }
f <- fit_of(X); M <- f$mi; C <- f$clr_scores; f$release(); rm(f); diag(M) <- 0

## ---- Task B diagnostics (label-free) ------------------------------------------------------
set.seed(20261103L)
Xp <- t(apply(X, 1, sample)); rownames(Xp) <- genes
fp <- fit_of(Xp); Mp <- fp$mi; Cp <- fp$clr_scores; fp$release(); rm(fp); diag(Mp) <- 0
R <- stats::cor(t(X)); R[is.na(R)] <- 0; diag(R) <- 0
Rp <- stats::cor(t(Xp)); Rp[is.na(Rp)] <- 0; diag(Rp) <- 0
dB <- rbind(cbind(measure = "mi_bspline", row_null_diag(M, Mp, "k562_truthgenes")),
            cbind(measure = "abs_r", row_null_diag(abs(R), abs(Rp), "k562_truthgenes")))
Ro <- R; diag(Ro) <- NA
dB$row_median_signed_r_med <- c(NA, stats::median(apply(Ro, 1, stats::median, na.rm = TRUE)))
dB$frac_signed_r_pos <- c(NA, mean(Ro > 0, na.rm = TRUE)); rm(Ro)
dB$n_profiles <- ncol(X)
utils::write.csv(dB, file.path(opt$out, "taskB_k562_diagnostics.csv"), row.names = FALSE)
rk <- rbind(rank_change(clr::clr_calibrate(abs(R), "normal", "stouffer"), C, "k562", ref = "clr_mi"))
rk$method <- "clr_absr"; utils::write.csv(rk, file.path(opt$out, "taskB_k562_rank_change.csv"), row.names = FALSE)
rm(R, Rp, Mp, Xp); invisible(gc(verbose = FALSE))
# label-free degree for the degree-scaled readout: permutation BH at q = 0.05 on the pooled null
ut <- upper.tri(C); obs <- C[ut]; nul <- sort(Cp[ut]); rm(Cp)
pv <- (1 + length(nul) - findInterval(obs, nul, left.open = TRUE)) / (1 + length(nul))   # P(null >= obs)
o <- order(pv); m <- length(pv); okk <- which(pv[o] <= seq_len(m) * 0.05 / m)
tau_bh <- if (length(okk)) min(obs[o][seq_len(max(okk))]) else Inf
n_bh <- if (is.finite(tau_bh)) sum(obs >= tau_bh & obs > 0) else 0L
# Degree for the degree-scaled readout (fixed before running; the BH set can be empty, see the smoke test):
# 1 + number of gene i's pairs in the global top 1% of CLR pairs (label-free).
q99 <- unname(stats::quantile(obs, 0.99)); E <- C >= q99 & C > 0; diag(E) <- FALSE; deg <- 1 + rowSums(E); rm(E, obs, nul, pv, o)
say(sprintf("permutation BH q=0.05 threshold on CLR %.3f (%d pairs pass); top-1%% CLR threshold %.3f; degree median %g, q90 %g, max %g",
            tau_bh, n_bh, q99, stats::median(deg), stats::quantile(deg, 0.9), max(deg)))
utils::write.csv(data.frame(tau_bh_q05 = tau_bh, n_pairs_pass_bh = n_bh, q99_clr = q99, deg_median = stats::median(deg),
  deg_q90 = unname(stats::quantile(deg, 0.9)), deg_max = max(deg)), file.path(opt$out, "taskB_k562_perm_bh.csv"), row.names = FALSE)

## ---- truths -------------------------------------------------------------------------------
RP_RE <- "^RP[LS][0-9]+[A-Z]?[0-9]*$|^RPLP[0-2]$|^RPSA$"       # analysis/clr_logpt_test.R
merge_sets <- function(s) {                                     # analysis/clr_logpt_test.R, verbatim
  if (length(s) == 0) return(s)
  s <- s[order(names(s))]
  key <- vapply(s, paste, "", collapse = " "); s <- s[!duplicated(key)]
  drop <- vapply(seq_along(s), function(i) any(vapply(seq_along(s)[-i], function(j)
    length(s[[j]]) > length(s[[i]]) && all(s[[i]] %in% s[[j]]), TRUE)), TRUE)
  s <- s[!drop]
  if (length(s) < 2) return(s)
  all_g <- sort(unique(unlist(s))); M <- vapply(s, function(g) all_g %in% g, logical(length(all_g)))
  inter <- crossprod(M); sz <- colSums(M); jac <- inter / (outer(sz, sz, "+") - inter)
  cl <- stats::cutree(stats::hclust(stats::as.dist(1 - jac), "complete"), h = 0.5)
  out <- lapply(split(seq_along(s), cl), function(ix) sort(unique(unlist(s[ix]))))
  names(out) <- vapply(split(names(s), cl), function(n) paste(n, collapse = "+"), "")
  out
}
none <- stats::setNames(character(0), character(0))
as_reg <- function(sets, cls) list(members = sets, class = stats::setNames(rep(cls, length(sets)), names(sets)))
cr <- corum_sets(genes); raw_m <- cr$members; raw_m <- raw_m[lengths(raw_m) >= 3L & lengths(raw_m) <= 200L]
mg <- merge_sets(lapply(raw_m, function(ix) genes[ix])); mg <- lapply(mg, function(g) sort(match(g, genes)))
rpfrac_mg <- vapply(mg, function(ix) mean(grepl(RP_RE, genes[ix])), 0)
sz <- lengths(raw_m)
TR <- list(CORUM_raw = as_reg(raw_m, "complex"), CORUM_merged = as_reg(mg, "complex"),
           CORUM_merged_noRP = as_reg(mg[rpfrac_mg < 0.3], "complex"),
           CORUM_small_3to6 = as_reg(raw_m[sz <= 6], "complex"), CORUM_large_ge20 = as_reg(raw_m[sz >= 20], "complex"),
           DoRothEA_AB = dorothea_sets(genes), CollecTRI = collectri_sets(genes))
bench <- lapply(TR, function(r) comembership_pairs(r, none, genes))
bench <- bench[vapply(bench, function(p) sum(p$label) > 0, NA)]
for (e in names(bench)) say(sprintf("truth %-18s %6d positives, %9d pairs, base rate %.4f", e, sum(bench[[e]]$label), length(bench[[e]]$label), mean(bench[[e]]$label)))
co <- matrix(FALSE, G, G); for (mm in raw_m) co[mm, mm] <- TRUE; diag(co) <- FALSE
big <- which(vapply(raw_m, function(mm) length(mm) * (length(mm) - 1) / 2 >= 10, NA))
rp_raw <- vapply(raw_m, function(ix) mean(grepl(RP_RE, genes[ix])), 0)

## ---- operators -----------------------------------------------------------------------------
selfloop <- function(A) { e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1; A }
lazy_P <- function(A, alpha = 0.5) Matrix::Diagonal(nrow(A), 1 - alpha) + alpha * Matrix::Matrix(A, sparse = TRUE)
perp <- function(A) { L <- A * log(A); L[!is.finite(L)] <- 0; exp(-rowSums(L)) }
TS <- c(1L, 2L, 4L, 9L, 16L)
OPS <- list(K10c50 = c(10, 50), K10c100 = c(10, 100), K10c200 = c(10, 200), K25c50 = c(25, 50), K40c50 = c(40, 50),
            K40c100 = c(40, 100), K80c100 = c(80, 100), K80c200 = c(80, 200), U50 = c(50, 50), U100 = c(100, 100))
sym <- function(Pt) { S <- (Pt + t(Pt)) / 2; diag(S) <- 0; S }
a3r <- function(Pt) { L <- log(Pt); L[!is.finite(L)] <- NA; diag(L) <- NA       # docs/plans/corum_analysis.md reference code
  mu <- rowMeans(L, na.rm = TRUE); n <- rowSums(!is.na(L))
  Z <- (L - mu) / sqrt(rowSums((L - mu)^2, na.rm = TRUE) / (n - 1)); Z[is.na(Z) | Z < 0] <- 0
  S <- (Z + t(Z)) / sqrt(2); diag(S) <- 0; S }
rowpct <- function(S0) { R <- S0; diag(R) <- NA
  R <- t(apply(R, 1, function(v) { k <- !is.na(v); o <- v; o[k] <- rank(v[k]) / sum(k); o }))
  S <- (R + t(R)) / 2; S[is.na(S)] <- 0; diag(S) <- 0; S }
minr <- function(Pt) { S <- pmin(Pt, t(Pt)); diag(S) <- 0; S }
degr <- function(Pt) { D <- deg * Pt; S <- (D + t(D)) / 2; diag(S) <- 0; S }

hub_stratify <- function(S, H, k = 10L) {                       # baselines.R, verbatim
  S[!is.finite(S)] <- 0; ut <- which(upper.tri(S)); h <- H[ut]
  g <- cut(h, unique(stats::quantile(h, seq(0, 1, length.out = k + 1L))), include.lowest = TRUE)
  r <- stats::ave(S[ut], g, FUN = function(x) rank(x) / length(x))
  R <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S)); R[ut] <- r; R + t(R)
}
# hub covariate of baselines.R: stationary distribution of the reference lazy walk (K10, cap 50), teleport 1e-3
P10 <- lazy_P(selfloop(clr:::.softmax_rows(C, keff = 10, cap = 50L)))
pi <- rep(1 / G, G); eps <- 1e-3
for (it in 1:20000) { q <- (1 - eps) * as.numeric(pi %*% P10) + eps / G; if (sum(abs(q - pi)) < 1e-13) break; pi <- q }
H <- outer(q, q, "+"); rm(P10)
HUBT <- c("CORUM_raw", "CORUM_merged", "DoRothEA_AB", "CollecTRI")
met <- list(); SC9 <- list(); MAC <- list(); PERP <- list(); HB <- list()
score_all <- function(name, S, t, keep9 = FALSE) {
  for (e in names(bench)) { b <- bench[[e]]; s <- S[cbind(b$i, b$j)]; s[!is.finite(s)] <- 0
    met[[length(met) + 1L]] <<- data.frame(variant = name, t = t, truth = e, aupr = ap_fast(s, b$label), base_rate = mean(b$label))
    if (keep9) SC9[[name]][[e]] <<- s }
  if (keep9) { x <- regulon_ap_table(S, list(members = raw_m[big]), co, none, genes); x$k <- big
    x$rp_dominated <- rp_raw[big] >= 0.3; x$variant <- name; MAC[[name]] <<- x
    Sh <- hub_stratify(S, H)
    for (e in intersect(HUBT, names(bench))) { b <- bench[[e]]
      HB[[length(HB) + 1L]] <<- data.frame(truth = e, variant = name, aupr_hubstrat = ap_fast(Sh[cbind(b$i, b$j)], b$label)) } }
}
score_all("CLR", C, 0L, TRUE); score_all("CLR_rowpct", rowpct(C), 0L, TRUE)
for (op in names(OPS)) {
  t0 <- Sys.time(); kc <- OPS[[op]]
  A <- selfloop(clr:::.softmax_rows(C, keff = kc[1], cap = as.integer(kc[2])))
  pa <- perp(A)
  PERP[[op]] <- data.frame(variant = op, keff = kc[1], cap = kc[2], perp_med = stats::median(pa),
                           perp_q10 = unname(stats::quantile(pa, 0.1)), frac_rows_uniform = mean(rowSums(A > 0) <= kc[1]))
  P <- lazy_P(A); rm(A)
  Pt <- NULL
  for (t in seq_len(max(TS))) {
    Pt <- if (is.null(Pt)) as.matrix(P) else as.matrix(Pt %*% P)
    if (!(t %in% TS)) next
    score_all(op, sym(Pt), t, keep9 = (t == 9L))
    if (op == "K10c50" && t == 9L) {
      score_all("K10c50_min", minr(Pt), t, TRUE); score_all("K10c50_rowpct", rowpct(Pt), t, TRUE)
      score_all("K10c50_degree", degr(Pt), t, TRUE); score_all("K10c50_a3r", a3r(Pt), t, TRUE)
    }
    if (op %in% c("K40c100", "U100") && t == 9L) score_all(paste0(op, "_a3r"), a3r(Pt), t, TRUE)
  }
  rm(P, Pt); invisible(gc(verbose = FALSE))
  say(sprintf("operator %-8s done (%.0f s)", op, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
utils::write.csv(do.call(rbind, HB), file.path(opt$out, "t9_hubstrat_EXPLORATORY.csv"), row.names = FALSE)
met <- do.call(rbind, met); met$enrichment <- met$aupr / met$base_rate
utils::write.csv(met, file.path(opt$out, "pooled_all_t.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, PERP), file.path(opt$out, "row_perplexity.csv"), row.names = FALSE)

## ---- t = 9 table with paired delete-half gene jackknife (seed 20261016, as baselines.R) ------------
set.seed(20261016L)
halves <- lapply(bench, function(b) lapply(seq_len(opt$jack), function(r) {
  keep <- sample(b$genes_in, floor(length(b$genes_in) / 2)); inh <- logical(G); inh[keep] <- TRUE
  which(inh[b$i] & inh[b$j]) }))
vars <- names(SC9); nc <- max(1L, min(opt$threads, parallel::detectCores() - 1L))
jk <- do.call(rbind, lapply(names(bench), function(e) { lab <- bench[[e]]$label
  J <- do.call(cbind, lapply(vars, function(v) unlist(parallel::mclapply(halves[[e]], function(ix) ap_fast(SC9[[v]][[e]][ix], lab[ix]), mc.cores = nc))))
  colnames(J) <- vars
  full <- vapply(vars, function(v) ap_fast(SC9[[v]][[e]], lab), 0)
  do.call(rbind, lapply(vars, function(v) { dK <- J[, v] - J[, "K10c50"]; dC <- J[, v] - J[, "CLR"]
    data.frame(truth = e, variant = v, aupr = full[[v]], enrichment = full[[v]] / mean(lab),
               d_vs_K10 = full[[v]] - full[["K10c50"]], sd_vs_K10 = stats::sd(dK),
               d_vs_CLR = full[[v]] - full[["CLR"]], sd_vs_CLR = stats::sd(dC)) })) }))
jk$lo_vs_K10 <- jk$d_vs_K10 - 1.96 * jk$sd_vs_K10; jk$hi_vs_K10 <- jk$d_vs_K10 + 1.96 * jk$sd_vs_K10
jk$lo_vs_CLR <- jk$d_vs_CLR - 1.96 * jk$sd_vs_CLR; jk$hi_vs_CLR <- jk$d_vs_CLR + 1.96 * jk$sd_vs_CLR
utils::write.csv(jk, file.path(opt$out, "t9_jackknife_EXPLORATORY.csv"), row.names = FALSE)

## ---- per-complex (macro) AP by size bin --------------------------------------------------------
mac <- do.call(rbind, MAC); mac$bin <- cut(mac$size, c(0, 9, 24, 49, Inf), labels = c("3-9", "10-24", "25-49", ">=50"))
utils::write.csv(mac, file.path(opt$out, "t9_macro_per_complex_EXPLORATORY.csv"), row.names = FALSE)
ref <- mac[mac$variant == "CLR", c("regulon", "ap")]; names(ref)[2] <- "ap_clr"
k10 <- mac[mac$variant == "K10c50", c("regulon", "ap")]; names(k10)[2] <- "ap_k10"
mm <- merge(merge(mac, ref), k10)
ms <- do.call(rbind, lapply(split(mm, list(mm$variant, mm$bin), drop = TRUE), function(d) data.frame(variant = d$variant[1], bin = d$bin[1],
  n = nrow(d), mean_ap = mean(d$ap, na.rm = TRUE), frac_beats_CLR = mean(d$ap > d$ap_clr, na.rm = TRUE),
  frac_beats_K10 = mean(d$ap > d$ap_k10, na.rm = TRUE), median_ratio_to_CLR = stats::median(d$ap / d$ap_clr, na.rm = TRUE))))
ms2 <- do.call(rbind, lapply(split(mm, mm$variant), function(d) data.frame(variant = d$variant[1], bin = "all", n = nrow(d),
  mean_ap = mean(d$ap, na.rm = TRUE), frac_beats_CLR = mean(d$ap > d$ap_clr, na.rm = TRUE), frac_beats_K10 = mean(d$ap > d$ap_k10, na.rm = TRUE),
  median_ratio_to_CLR = stats::median(d$ap / d$ap_clr, na.rm = TRUE))))
utils::write.csv(rbind(ms2, ms), file.path(opt$out, "t9_macro_summary_EXPLORATORY.csv"), row.names = FALSE)

chk <- jk[jk$truth == "CORUM_raw" & jk$variant %in% c("CLR", "K10c50", "K10c50_a3r"), c("variant", "aupr", "enrichment")]
chk$documented <- c(CLR = 0.0847135108446739, K10c50 = 0.0465396136302405, K10c50_a3r = 0.105435082561767)[chk$variant]
utils::write.csv(chk, file.path(opt$out, "repro_check.csv"), row.names = FALSE); print(chk)
say("done")
