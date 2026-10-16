#!/usr/bin/env Rscript
# Pilot for the label-free multiscale algorithm (paper plan, theme 4): are there
# plateaus in the spectrum of the attention walk, and does keff decide which
# levels are visible? No bootstrap yet.
#
# For each dataset and keff in {3, 5, 10}: softmax operator A (cap 50), lazy walk
# P = (1 - a) I + a A, stationary distribution pi (teleportation eps = 1e-3), and
# its reversibilization P_rev = (P + Pi^-1 P' Pi) / 2, which is similar to a
# symmetric matrix. One eigendecomposition then gives, for any real t:
#   N(t)   = trace(P_rev^t) = sum_k mu_k^t   effective number of modules
#   r_i(t) = (P_rev^t)_ii                    return probability of gene i
#   k_i(t) = G pi_i / r_i(t)                 effective module size of gene i (genes)
# Regimes are plateaus of log N against log t (local slope near 0); plateau
# height is the module count. Null: the same operator built on CLR scores
# permuted among gene pairs (keff = 10).
# Hub check: the degree-only score pi_i + pi_j scored against the truths, to see
# how much of the large-t AUPR is convergence to hubs.
#
# Usage (repo root):  Rscript analysis/multiscale_pilot.R [--out results/multiscale_pilot]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(out = "results/multiscale_pilot", alpha = 0.5, eps = 1e-3, seed = 20261015L, p_max = 1e-3,
            h5 = "data/perturbseq/K562_essential_normalized_bulk_01.h5ad", keffs = "3,5,10")
i <- 1L
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1L]; i <- i + 2L }
keffs <- as.numeric(strsplit(opt$keffs, ",")[[1]])
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
logf <- file.path(opt$out, "multiscale_pilot.log")
say <- function(...) { m <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...))
                       cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }
suppressPackageStartupMessages({ library(clr); library(Matrix) })
source("analysis/regulondb.R"); source("analysis/regulons.R"); source("analysis/human_truths.R")
say("options: ", paste(names(opt), unlist(opt), sep = "=", collapse = " "))

tgrid <- exp(seq(0, log(4096), length.out = 121))

operator <- function(C, keff) {
  A <- clr:::.softmax_rows(C, keff = keff, cap = 50L); e <- which(rowSums(A) == 0); A[cbind(e, e)] <- 1
  (1 - opt$alpha) * diag(nrow(A)) + opt$alpha * A
}
stationary <- function(P) {
  G <- nrow(P); p <- rep(1 / G, G); Ps <- Matrix::Matrix(P, sparse = TRUE)
  for (it in 1:20000) { q <- (1 - opt$eps) * as.numeric(p %*% Ps) + opt$eps / G
                        if (sum(abs(q - p)) < 1e-13) break; p <- q }
  q / sum(q)
}
spectrum <- function(P, label) {
  G <- nrow(P); pi <- stationary(P)
  Pe <- (1 - opt$eps) * P + opt$eps / G
  M <- sqrt(pi) * Pe * rep(1 / sqrt(pi), each = G)          # Pi^1/2 Pe Pi^-1/2
  M <- (M + t(M)) / 2; rm(Pe)                              # = Pi^1/2 P_rev Pi^-1/2
  ev <- eigen(M, symmetric = TRUE); rm(M)
  if (min(ev$values) < -1e-8) say(sprintf("%s: %d negative eigenvalues clipped (min %.2g)", label,
                                         sum(ev$values < -1e-8), min(ev$values)))
  mu <- pmax(ev$values, 0); U2 <- ev$vectors^2; rm(ev)
  powt <- outer(mu, tgrid, `^`)                            # G x T
  N <- colSums(powt); R <- U2 %*% powt                     # return probabilities, G x T
  logk <- log(G * pi) - log(R)
  slope <- c(NA, diff(log(N)) / diff(log(tgrid)))
  q <- apply(logk, 2, stats::quantile, c(0.25, 0.5, 0.75))
  list(pi = pi, curve = data.frame(op = label, t = tgrid, N = N, slope = slope,
    k_q25 = exp(q[1, ]), k_med = exp(q[2, ]), k_q75 = exp(q[3, ]),
    k_med_slope = c(NA, diff(q[2, ]) / diff(log(tgrid)))))
}
plateaus <- function(cur, thr = 0.1) {                     # runs of |slope| < thr, N >= 2
  ok <- !is.na(cur$slope) & abs(cur$slope) < thr & cur$N >= 2
  r <- rle(ok); end <- cumsum(r$lengths); start <- end - r$lengths + 1
  k <- which(r$values & r$lengths >= 3)
  if (!length(k)) return(NULL)
  data.frame(op = cur$op[1], t_from = cur$t[start[k]], t_to = cur$t[end[k]],
             log_t_width = log(cur$t[end[k]] / cur$t[start[k]]),
             N_mid = cur$N[(start[k] + end[k]) %/% 2])
}
hub_check <- function(pi, bench, label) do.call(rbind, lapply(names(bench), function(nm) {
  p <- bench[[nm]]
  data.frame(data = label, truth = nm, base_rate = mean(p$label),
             aupr_hub_sum = ap_fast(pi[p$i] + pi[p$j], p$label),
             aupr_hub_prod = ap_fast(pi[p$i] * pi[p$j], p$label)) }))
null_scores <- function(C) {                               # CLR values permuted among pairs
  set.seed(opt$seed); ut <- which(upper.tri(C)); N <- C; N[ut] <- sample(C[ut])
  N[lower.tri(N)] <- t(N)[lower.tri(N)]; N
}
run_dataset <- function(C, bench, label) {
  curves <- list(); pls <- list(); hub <- NULL
  for (kf in keffs) {
    P <- operator(C, kf); s <- spectrum(P, sprintf("%s keff=%g", label, kf)); rm(P)
    curves[[length(curves) + 1L]] <- s$curve; pls[[length(pls) + 1L]] <- plateaus(s$curve)
    if (kf == 10) hub <- hub_check(s$pi, bench, label)
    say(sprintf("%s keff=%g: N(1) = %.0f, N(64) = %.1f, N(1024) = %.2f; median k(8) = %.1f, k(64) = %.1f genes",
      label, kf, s$curve$N[1], approx(tgrid, s$curve$N, 64)$y, approx(tgrid, s$curve$N, 1024)$y,
      approx(tgrid, s$curve$k_med, 8)$y, approx(tgrid, s$curve$k_med, 64)$y))
  }
  P <- operator(null_scores(C), 10); s <- spectrum(P, sprintf("%s keff=10 permuted-CLR null", label)); rm(P)
  curves[[length(curves) + 1L]] <- s$curve; pls[[length(pls) + 1L]] <- plateaus(s$curve)
  cur <- do.call(rbind, curves); cur$data <- label
  list(curves = cur, plateaus = do.call(rbind, pls), hub = hub)
}

## ---- E. coli (as granularity.R) ------------------------------------------------------------
m3d <- "data/E_coli_v4_Build_6"; base <- basename(m3d)
raw <- utils::read.delim(file.path(m3d, paste0("avg_", base, "_exps466probes4297.tab")), check.names = FALSE)
X <- as.matrix(raw[, -1]); rownames(X) <- vapply(strsplit(raw[[1]], "_"), function(p) p[length(p) - 1L], "")
X <- X[apply(X, 1, function(x) diff(range(x))) > 0, ]; genes <- rownames(X); G <- length(genes)
to_bn <- make_bnumber_mapper("analysis/ecoli_k12_genes.tsv")
rtab <- build_regulon_table("data/RegulonDBExtract", to_bn)
operon_of <- build_operon_map("data/RegulonDBExtract", to_bn)
reg_all <- make_regulons(rtab, genes, conf_keep = c("C", "S", "W", "?"), max_size = 2000)
pick <- function(keep) list(members = reg_all$members[keep], class = reg_all$class[keep])
tf_reg <- pick(reg_all$class == "TF" & lengths(reg_all$members) <= 500 & lengths(reg_all$members) >= 5)
sig_reg <- pick(reg_all$class == "sigma" & lengths(reg_all$members) >= 5 &
                  !names(reg_all$members) %in% c("sigma:sigma70", "sigma:sigma38"))
op <- operon_of[names(operon_of) %in% genes]
op_sets <- lapply(split(match(names(op), genes), op), unique); op_sets <- op_sets[lengths(op_sets) >= 2]
op_reg <- list(members = op_sets, class = stats::setNames(rep("operon", length(op_sets)), names(op_sets)))
none <- stats::setNames(character(0), character(0))
bench_e <- list(`1 operons` = comembership_pairs(op_reg, none, genes),
                `2 TF regulons` = comembership_pairs(tf_reg, operon_of, genes),
                `3 sigmulons` = comembership_pairs(sig_reg, operon_of, genes))
f <- ClrAttention$new(X)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")
C <- f$clr_scores; f$release(); rm(f, X); invisible(gc())
say("E. coli M3D: ", G, " genes")
res_e <- run_dataset(C, bench_e, "E. coli M3D"); rm(C); invisible(gc())

## ---- K562 essential compendium, network half A (as granularity.R) ---------------------------
set.seed(opt$seed)
rd <- function(p) rhdf5::h5read(opt$h5, p)
cat_of <- function(grp, field) { x <- rd(paste0(grp, "/", field))
  if (is.list(x)) return(as.character(x$categories)[x$codes + 1L])
  cats <- tryCatch(rd(paste0(grp, "/__categories/", field)), error = function(e) NULL)
  if (is.null(cats)) as.character(x) else as.character(cats)[x + 1L] }
Z <- rd("X"); hg <- cat_of("/var", "gene_name"); rownames(Z) <- hg
pert <- cat_of("/obs", "gene_transcript"); colnames(Z) <- pert
pval <- as.numeric(rd("/obs/energy_test_p_value"))
K <- pert[pval < opt$p_max & !grepl("non-targeting", pert)]; K <- K[sample(length(K))]
A <- sort(K[seq_len(floor(length(K) / 2))]); B <- setdiff(K, A)
Xh <- Z[, A]; rm(Z); Xh[!is.finite(Xh)] <- 0
Xh <- Xh[apply(Xh, 1, function(x) diff(range(x))) > 0, ]; hgenes <- rownames(Xh)
pairs_h <- function(r) comembership_pairs(r, none, hgenes)
bench_h <- list(`1 CORUM complexes` = pairs_h(corum_sets(hgenes)),
                `2 KD-response modules` = pairs_h(response_sets(hgenes, path = opt$h5, use = B)),
                `3 DoRothEA A-B` = pairs_h(dorothea_sets(hgenes)),
                `4 DoRothEA A-C` = pairs_h(dorothea_sets(hgenes, c("A", "B", "C"))),
                `5 CollecTRI` = pairs_h(collectri_sets(hgenes)))
fh <- ClrAttention$new(Xh)$estimate_mi(bins = "hg", transform = "none")$calibrate(method = "normal", combine = "stouffer")
C <- fh$clr_scores; fh$release(); rm(fh, Xh); invisible(gc())
say("K562 essential compendium: ", nrow(C), " genes")
res_h <- run_dataset(C, bench_h, "K562 essential compendium"); rm(C); invisible(gc())

## ---- output -----------------------------------------------------------------------------------
utils::write.csv(rbind(res_e$curves, res_h$curves), file.path(opt$out, "spectrum_curves.csv"), row.names = FALSE)
pl <- rbind(res_e$plateaus, res_h$plateaus)
if (!is.null(pl)) utils::write.csv(pl, file.path(opt$out, "plateaus.csv"), row.names = FALSE)
hub <- rbind(res_e$hub, res_h$hub); hub$hub_ratio <- hub$aupr_hub_sum / hub$base_rate
utils::write.csv(hub, file.path(opt$out, "hub_baseline.csv"), row.names = FALSE)
say("plateaus (|d log N / d log t| < 0.1 over >= 3 grid points, N >= 2):")
if (is.null(pl)) say("  none") else for (k in seq_len(nrow(pl)))
  say(sprintf("  %-45s t = %7.1f to %7.1f  (width %.2f in log t)  N = %.1f", pl$op[k], pl$t_from[k], pl$t_to[k],
              pl$log_t_width[k], pl$N_mid[k]))
say("hub baseline (pi_i + pi_j, keff = 10 operator):")
for (k in seq_len(nrow(hub))) say(sprintf("  %-26s %-22s AUPR %.4f (x%.2f base rate; product %.4f)",
  hub$data[k], hub$truth[k], hub$aupr_hub_sum[k], hub$hub_ratio[k], hub$aupr_hub_prod[k]))
say("done")
