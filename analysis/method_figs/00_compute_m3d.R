# Worked-example numbers for the eight method figures (docs/checkpoint_collins_src/figures/method_*.png). Everything the figures print is computed here from the real
# M3D matrix (4297 genes x 466 profiles, the matrix of the tcaCLR / tcaCLR2 evaluation) with the project's own functions:
#   fit (B-spline MI, row-z CLR score)          analysis/tcaclr_depth/grid_lib.R gl_fit via adaptive_lib.R ad_fit
#   operator P = 0.5 I + 0.5 A (keff 10, cap 50) gl_make_P
#   walk P^9, second CLR                         gl_power_seq, gl_zmat
#   permutation edge (10 per-gene permutations)  ad_permute_genes, ad_Psym, ad_eigs (the seeds of results/tcaclr2/edge_*.tsv)
#   global temperature, 13-temperature ladder    analysis/global_temp_ladder/gt_lib.R
# No truth (regulon, complex) is read. Output: analysis/method_figs/cache/m3d_example.rds (small objects only).
# Run from the repo root:  nice -n 19 Rscript analysis/method_figs/00_compute_m3d.R
suppressPackageStartupMessages(source("analysis/global_temp_ladder/gt_lib.R"))
dir.create("analysis/method_figs/cache", showWarnings = FALSE, recursive = TRUE)
THREADS <- 8L; KEFF <- 10; CAP <- 50L; T9 <- 9L
say <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

X <- gl_truth_matrix("m3d_tf_regulons"); G <- nrow(X); NS <- ncol(X); genes <- rownames(X)
gtab <- utils::read.delim("analysis/ecoli_k12_genes.tsv", stringsAsFactors = FALSE); syms <- gtab$symbol[match(genes, gtab$bnumber)]; syms[is.na(syms) | !nzchar(syms)] <- genes[is.na(syms) | !nzchar(syms)]
say(sprintf("M3D: %d genes x %d profiles; duplicated gene names: %d", G, NS, sum(duplicated(genes))))
f <- ad_fit(X, THREADS); S <- f$clr_scores; D <- f$mi; diag(D) <- 0
say("fit done")

## ---- CLR by hand (same arithmetic as clr_core.cpp) and check against the package score ----------------------------------------------
mu <- rowSums(D) / G; Dc <- D - mu; sdv <- sqrt(rowSums(Dc^2) / (G - 1)); Zraw <- Dc / sdv; rm(Dc)   # unclipped row z (diagonal zeroed MI enters mu and sd)
Zc <- pmax(Zraw, 0); S_chk <- (Zc + t(Zc)) / sqrt(2); diag(S_chk) <- 0
stopifnot(max(abs(S_chk - S)) < 1e-9); rm(S_chk, Zc); say("hand CLR equals package CLR")

## ---- candidate lists, per-gene temperature (the bisection of clr:::.softmax_rows), effective partners ---------------------------------
Tm <- gt_topmat(S, CAP)
ncand <- rowSums(!is.na(Tm$idx))
neff_of <- function(s, tau) { w <- exp((s - max(s)) / tau); w <- w / sum(w); exp(-sum(w[w > 0] * log(w[w > 0]))) }
tau_row <- function(s, keff = KEFF) { lo <- log(1e-4); hi <- log(1e4); for (it in 1:60) { mid <- (lo + hi) / 2; if (neff_of(s, exp(mid)) < keff) lo <- mid else hi <- mid }; exp((lo + hi) / 2) }
tau_loc <- rep(NA_real_, G); eff_loc <- rep(NA_real_, G)
for (g in seq_len(G)) { s <- Tm$sc[g, !is.na(Tm$sc[g, ])]; if (length(s) == 0L) { eff_loc[g] <- 1; next }
  if (length(s) <= KEFF) { eff_loc[g] <- length(s); next }   # at most keff candidates: uniform weights (as .softmax_rows)
  tau_loc[g] <- tau_row(s); eff_loc[g] <- neff_of(s, tau_loc[g]) }
A <- gl_operator_A(S, KEFF, CAP)     # the package operator
set.seed(1); chk <- sample(which(!is.na(tau_loc)), 25, replace = FALSE)
for (g in chk) { s <- Tm$sc[g, !is.na(Tm$sc[g, ])]; w <- exp((s - max(s)) / tau_loc[g]); w <- w / sum(w); stopifnot(max(abs(A[g, Tm$idx[g, seq_along(s)]] - w)) < 1e-9) }
say(sprintf("per-gene tau check on 25 rows against the package operator: ok. rows with > %d candidates: %d of %d (tau solved); <= %d candidates (uniform): %d", KEFF, sum(!is.na(tau_loc)), G, KEFF, sum(is.na(tau_loc))))

## ---- walk, second CLR --------------------------------------------------------------------------------------------------------------
P <- gl_make_P(S, KEFF, CAP)
stopifnot(max(abs(as.matrix(P) - (0.5 * diag(G) + 0.5 * A))) < 1e-12)
P9 <- gl_power_seq(P, T9); say("P^9 done")
logP <- log(P9); logP[!is.finite(logP)] <- NA_real_; diag(logP) <- NA_real_
mu2 <- rowMeans(logP, na.rm = TRUE); n2 <- rowSums(!is.na(logP)); sd2 <- sqrt(rowSums((logP - mu2)^2, na.rm = TRUE) / (n2 - 1))
Z2 <- gl_zmat(P9)                      # row z over reachable off-diagonal entries; unreachable and diagonal = 0
chk2 <- (logP - mu2) / sd2; chk2[is.na(chk2)] <- 0; stopifnot(max(abs(chk2 - Z2)) < 1e-9); rm(chk2)
S2 <- (pmax(Z2, 0) + pmax(t(Z2), 0)) / sqrt(2); diag(S2) <- 0
aclr <- (P9 + t(P9)) / 2; diag(aclr) <- 0

## ---- neighbourhood: gene i = lexA, its five highest CLR partners, three genes ranked highest by the second CLR that are not among the five, two background genes ----
i <- match("b4043", genes)   # lexA; stopifnot(!is.na(i))
r1 <- S[i, ]; r1[i] <- -Inf; top5 <- order(r1, decreasing = TRUE)[1:5]
r2 <- S2[i, ]; r2[c(i, top5)] <- -Inf; walk3 <- order(r2, decreasing = TRUE)[1:3]
# two background genes: drawn at random (fixed seed) among the genes that have no positive CLR score with lexA, so the clipping step has something to clip
set.seed(7); bg2 <- sample(setdiff(which(S[i, ] == 0), c(i, top5, walk3)), 2)
Nb <- c(i, top5, walk3, bg2); nm <- syms[Nb]
say("neighbourhood: ", paste(nm, collapse = " "))

## ---- genes for the local-temperature picture: lexA and three genes at the 10th/50th/90th percentile of tau among genes with a full 50 candidates ------------
full50 <- which(!is.na(tau_loc) & ncand == CAP)
qs <- stats::quantile(tau_loc[full50], c(0.1, 0.5, 0.9), type = 1)
ex_g <- vapply(qs, function(q) full50[which.min(abs(tau_loc[full50] - q))], 1L)
ex_genes <- unique(c(i, ex_g)); say("local-tau genes: ", paste(syms[ex_genes], sprintf("tau=%.3f", tau_loc[ex_genes]), collapse = "; "))
tau_v2 <- gt_tau_for_eff(Tm)          # one global temperature, mean effective partners 10
w_g <- gt_w(Tm, tau_v2); has <- rowSums(w_g) > 0; wn <- w_g / pmax(rowSums(w_g), 1e-300)
eff_glob <- rep(NA_real_, G); eff_glob[has] <- gt_effv(Tm, tau_v2)
row_info <- lapply(ex_genes, function(g) { k <- sum(!is.na(Tm$sc[g, ])); s <- Tm$sc[g, seq_len(k)]
  wl <- exp((s - max(s)) / tau_loc[g]); wl <- wl / sum(wl); wg <- exp((s - max(s)) / tau_v2); wg <- wg / sum(wg)
  list(gene = syms[g], idx = Tm$idx[g, seq_len(k)], score = s, tau = tau_loc[g], w_loc = wl, w_glob = wg, eff_loc = neff_of(s, tau_loc[g]), eff_glob = neff_of(s, tau_v2),
       spread = c(max = max(s), q10th = s[10], q50th = s[min(50, k)], sd_all = stats::sd(S[g, ]))) })
names(row_info) <- syms[ex_genes]
say(sprintf("tau_v2 = %.4f (stored M3D table: see below); global effective partners: mean %.3f median %.2f q10 %.2f q90 %.2f min %.2f max %.2f",
            tau_v2, mean(eff_glob, na.rm = TRUE), median(eff_glob, na.rm = TRUE), quantile(eff_glob, .1, na.rm = TRUE), quantile(eff_glob, .9, na.rm = TRUE), min(eff_glob, na.rm = TRUE), max(eff_glob, na.rm = TRUE)))
tt <- utils::read.csv("results/global_temp_ladder/m3d/tau_table.csv"); stored_v2 <- tt$tau[tt$op == "V2"]
stopifnot(abs(tau_v2 - stored_v2) < 1e-6 * stored_v2); say(sprintf("tau_v2 equals the stored table value %.6f", stored_v2))
ect <- utils::read.csv("results/global_temp_ladder/ecomac/tau_table.csv"); ec_v2 <- ect[ect$op == "V2", ]

## ---- permutations and the edge ----------------------------------------------------------------------------------------------------
KV <- 40L
ev_obs <- ad_eigs(ad_Psym(P), KV)$values
perm_ev <- matrix(NA_real_, KV, ED$perms)
for (p in seq_len(ED$perms)) { fp <- ad_fit(ad_permute_genes(X, ED$seed + 7000L + p), THREADS)
  perm_ev[, p] <- ad_eigs(ad_Psym(gl_make_P(fp$clr_scores, KEFF, CAP)), KV)$values; fp$release(); rm(fp); invisible(gc(FALSE)); say("permutation ", p) }
edge_perm <- perm_ev[2, ]; edge <- mean(edge_perm); tstar <- log(ED$eps) / log(edge); t2 <- ed_t2(edge)
st <- utils::read.delim("results/tcaclr2/edge_m3d_tf_regulons.tsv")
say(sprintf("edge %.6f (stored %.6f), range %.6f-%.6f (stored %.6f-%.6f), t2 %d (stored %d), trivial eigenvalue min over perms %.10f",
            edge, st$edge, min(edge_perm), max(edge_perm), st$edge_min, st$edge_max, t2, st$t2, min(perm_ev[1, ])))
stopifnot(abs(edge - st$edge) < 1e-6, t2 == st$t2)
edge_files <- c(m3d = "results/tcaclr2/edge_m3d_tf_regulons.tsv", ecomac = "results/tcaclr2_scoreboard/edge_ecomac.tsv", k562 = "results/tcaclr2/edge_k562_corum.tsv",
                hepg2 = "results/tcaclr2/edge_armA_cp_human.tsv", kemmeren = "results/tcaclr2/edge_armB_cyc2008.tsv", yeast = "results/tcaclr2/edge_net4_cp_yeast.tsv")
edge_tab <- do.call(rbind, lapply(names(edge_files), function(d) { x <- utils::read.delim(edge_files[[d]]); data.frame(dataset = d, G = x$G, N = x$N, edge = x$edge, edge_min = x$edge_min, edge_max = x$edge_max, t2 = x$t2, tstar = x$tstar, perms = x$perms, eps = x$eps) }))
print(edge_tab)

## ---- multi-scale: 13 temperatures x second CLR at depth 9, percentile / z of every pair within the temperature ---------------------------
ut <- upper.tri(S); Npairs <- sum(ut)
# position of each neighbourhood pair in the column-major upper triangle (the order of S[ut])
pos <- matrix(0L, length(Nb), length(Nb)); cum <- c(0, cumsum(seq_len(G - 1L)))   # column c holds pairs (1..c-1, c): offset cum[c-1]
for (a in seq_along(Nb)) for (b in seq_along(Nb)) if (a != b) { r <- min(Nb[a], Nb[b]); cc <- max(Nb[a], Nb[b]); pos[a, b] <- cum[cc - 1L] + r }
Sut <- S[ut]; for (a in seq_along(Nb)) for (b in seq_along(Nb)) if (a != b) stopifnot(Sut[pos[a, b]] == S[Nb[a], Nb[b]]); rm(Sut)
pct <- array(NA_real_, c(length(Nb), length(Nb), length(GTL$taus))); zsc <- pct; sc_tau <- pct; eff_tau <- numeric(length(GTL$taus))
for (k in seq_along(GTL$taus)) { tau <- GTL$taus[k]
  Pk <- as.matrix(gt_make_P(Tm, tau)); Pk <- gl_power_seq(Pk, T9); Zk <- gl_zmat(Pk); rm(Pk)
  Sk <- (pmax(Zk, 0) + pmax(t(Zk), 0)) / sqrt(2); rm(Zk); v <- Sk[ut]
  rk <- rank(v, ties.method = "average") / Npairs; zk <- (v - mean(v)) / stats::sd(v)
  for (a in seq_along(Nb)) for (b in seq_along(Nb)) if (a != b) { pct[a, b, k] <- rk[pos[a, b]]; zsc[a, b, k] <- zk[pos[a, b]]; sc_tau[a, b, k] <- v[pos[a, b]] }
  eff_tau[k] <- mean(gt_effv(Tm, tau)); rm(Sk, v, rk, zk); invisible(gc(FALSE)); say(sprintf("ladder log2 tau %+.1f: mean effective partners %.2f", log2(tau), eff_tau[k])) }
ms1 <- apply(pct, c(1, 2), mean); ms3_idx <- which(round(log2(GTL$taus), 2) %in% round(GTL$ms3_log2, 2)); ms1_8 <- apply(pct[, , ms3_idx], c(1, 2), mean); ms2 <- apply(zsc, c(1, 2), mean)
stopifnot(length(ms3_idx) == 8L)

## ---- save small objects ------------------------------------------------------------------------------------------------------------
ex <- list(dataset = "M3D", G = G, N = NS, names = nm, Nb = Nb, i = i, top5 = top5, walk3 = walk3, bg2 = bg2,
  D_sub = D[Nb, Nb], mi_row_i = D[i, ], mu_i = mu[i], sd_i = sdv[i], mu_all = mu[Nb], sd_all = sdv[Nb], Zraw_sub = Zraw[Nb, Nb], S_sub = S[Nb, Nb], S_row_i = S[i, ],
  tau_loc = tau_loc, eff_loc = eff_loc, ncand = ncand, tau_v2 = tau_v2, eff_glob = eff_glob, row_info = row_info, ex_genes = syms[ex_genes],
  A_sub = A[Nb, Nb], P9_sub = P9[Nb, Nb], aclr_sub = aclr[Nb, Nb], logP_sub = logP[Nb, Nb], mu2 = mu2[Nb], sd2 = sd2[Nb], n2 = n2[Nb], Z2_sub = Z2[Nb, Nb], S2_sub = S2[Nb, Nb],
  logp_row_i = logP[i, ], reach_i = n2[i], reach_median = median(n2), reach_min = min(n2),
  A_row_i_nnz = sum(A[i, ] > 0),
  ev_obs = ev_obs, perm_ev = perm_ev, edge = edge, edge_perm = edge_perm, tstar = tstar, t2 = t2, edge_tab = edge_tab,
  ladder_log2 = log2(GTL$taus), pct = pct, zsc = zsc, sc_tau = sc_tau, ms1 = ms1, ms1_8 = ms1_8, ms2 = ms2, ms3_log2 = GTL$ms3_log2, eff_tau = eff_tau, Npairs = Npairs,
  eco_v2 = ec_v2)
saveRDS(ex, "analysis/method_figs/cache/m3d_example.rds"); say("saved")
